// hopkinsheat.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
// 
//////////////////////////////////////////////////////////////////////

#include "hopkinsheat.h"
#include "jacrow.h"
#include "cell.h"
#include "laminarduct.h"
#include "geom.h"

#define Ns (v.gc->Ns())

/* For obtaining time-averaged heat transfer, we have to find the
   time-average of the following parameter:
 
   i * omg * f_x
   ------------
   (1 - f_x )


   Heat transfer is related to temperature difference, according to the
   Hopkins papers with:

   Qin = Sf*H*(Tw-T) - Q*dTwdx*m


*/

namespace tube{
  using tasystem::JacRow;
  using tasystem::JacCol;
  using rottfuncs::RottFuncs;
  using tasystem::var;

  namespace H{
    typedef  double d;
    namespace{ 
      d zeroheat_vert(d kappa,d rh){
        return 3*kappa/pow(rh,2);
      }
      d zeroheat_circ(d kappa,d rh){
        return 2*kappa/pow(rh,2);
      }
      d zeroheat_inviscid(d,d){
        return 0;
      }
    }
  } // namespace H
  

  HopkinsHeatSource::HopkinsHeatSource(const LaminarDuct& t):
    cshape(t.geom().shape())
  {
    TRACE(10,"HopkinsHeatSource::HopKinsHeatSource(Tube)");
    this->t=&t;			// Save pointer to tube
    setZeroFreq(cshape);
    if(t.geom().isBlApprox())
      rf=RottFuncs("blapprox");
    else
      rf=RottFuncs(t.geom().shape()); // Reinitialize thermoviscous functions with right shape
    TRACE(11,"Exiting redefinition of Rottfuncs");
  }
  void HopkinsHeatSource::setZeroFreq(const string& shape){
    TRACE(10,"HopkinsHeatSource::setZeroFreq()");
    if(shape.compare("vert")==0){
      zeroheatH_funptr=&H::zeroheat_vert;
      zeroheatQ=1.0/5.0;      
    }
    else if(shape.compare("circ")==0){
      zeroheatH_funptr=&H::zeroheat_circ;
      zeroheatQ=1.0/3.0;      
    }
    else if(shape.compare("vert")==0){
      zeroheatH_funptr=&H::zeroheat_vert;
      zeroheatQ=0;
    }
    else if(shape.compare("inviscid")==0){
      zeroheatH_funptr=&H::zeroheat_inviscid;
      zeroheatQ=0;
    }
    else
      {
        WARN("Warning: tube.geom.shape unknown for ZeroHeatH. Aborting...");
        abort();
      }
  }
  std::tuple<d,d,d> wfacs(d xi,d xj,d xk) {
    d dj=xj-xi;
    d dk=xk-xi;
    d denom=dk*pow(dj,2)-dj*pow(dk,2);    
    d Wi=(pow(dk,2)-pow(dj,2))/denom;
    d Wj=-pow(dk,2)/denom;
    d Wk=pow(dj,2)/denom;
    return std::make_tuple(Wi,Wj,Wk);
  }
  d HopkinsHeatSource::dTwdx(const Cell& v) const {
    TRACE(15,"HopkinsHeatSource::dTwdx()");
    d xi,xj,xk,Wi,Wj,Wk;
    // VARTRACE(45,v.Tw()(0));
    if(v.left()&&v.right()){
      xi=v.vx;
      xj=v.left()->vx;
      xk=v.right()->vx;
      std::tie(Wi,Wj,Wk)=wfacs(xi,xj,xk);
      return Wi*v.Tw()(0)+Wj*v.TwL()(0)+Wk*v.TwR()(0);    
    }
    else if(!v.left()&&v.right()){
      xi=v.vx;
      xj=v.right()->vx;
      xk=v.right()->right()->vx;
      std::tie(Wi,Wj,Wk)=wfacs(xi,xj,xk);
      return Wi*v.Tw()(0)+Wj*v.TwR()(0)+Wk*v.right()->TwR()(0);    
    }
    else{
      xi=v.vx;
      xj=v.left()->vx;
      xk=v.left()->left()->vx;
      std::tie(Wi,Wj,Wk)=wfacs(xi,xj,xk);
      return Wi*v.Tw()(0)+Wj*v.TwL()(0)+Wk*v.left()->TwL()(0);    
    }
  }
  vd HopkinsHeatSource::Qsf(const Cell& v) const{
    TRACE(5,"HopkinsHeatSource::heat(v)");
    
    dmat H=var(*v.gc,HeatTransferCoefH(v)).freqMultiplyMat();
    vd heat=-v.vSf*(H*v.T()());
    heat(0)+=v.vSf*(H(0,0)*v.Tw()(0)); // Time-averaged component subtracted

    dmat Q=var(*v.gc,HeatTransferCoefQ(v)).freqMultiplyMat();
    heat-=Q*(dTwdx(v)*0.5*(v.ml()()+v.mr()()) );   
    return heat;    
  }
  JacRow HopkinsHeatSource::dQsf(const Cell& v) const{
    TRACE(15,"HopkinsHeatSource::dQsf()");
    JacRow dQsf(-1,4);
    dmat H=var(*v.gc,HeatTransferCoefH(v)).freqMultiplyMat();
    dmat Q=var(*v.gc,HeatTransferCoefQ(v)).freqMultiplyMat();    

    dmat H11=zeros(Ns,Ns);
    H11(0,0)=H(0,0);
    // Obtain dTwdx
    d dTwdx=this->dTwdx(v);
    dQsf+=JacCol(v.T(),-v.vSf*H);

    if(t->hasSolid()){

      dQsf+=JacCol(v.Tw(),v.vSf*H11);

      vd m=0.5*(v.ml()()+v.mr()());
      // Only zero in right top
      dmat row(1,Ns,fillwith::zeros); row(0,0)=1;
      d xi,xj,xk,Wi,Wj,Wk;
      if(v.left()&&v.right()){
	xi=v.vx;
	xj=v.left()->vx;
	xk=v.right()->vx;
	std::tie(Wi,Wj,Wk)=wfacs(xi,xj,xk);
	dQsf+=JacCol(v.Tw(),-Wi*(Q*m)*row);
	dQsf+=JacCol(v.TwL(),-Wj*(Q*m)*row);
	dQsf+=JacCol(v.TwR(),-Wk*(Q*m)*row);      
      }
      else if(!v.left()&&v.right()){
	// Leftmost cell
	xi=v.vx;
	xj=v.right()->vx;
	xk=v.right()->right()->vx;
	std::tie(Wi,Wj,Wk)=wfacs(xi,xj,xk);
	dQsf+=JacCol(v.Tw(),-Wi*(Q*m)*row);
	dQsf+=JacCol(v.TwR(),-Wj*(Q*m)*row);
	dQsf+=JacCol(v.right()->TwR(),-Wk*(Q*m)*row);      
      }
      else{
	// Rightmost cell
	xi=v.vx;
	xj=v.left()->vx;
	xk=v.left()->left()->vx;
	std::tie(Wi,Wj,Wk)=wfacs(xi,xj,xk);
	dQsf+=JacCol(v.Tw(),-Wi*(Q*m)*row);
	dQsf+=JacCol(v.TwL(),-Wj*(Q*m)*row);
	dQsf+=JacCol(v.left()->TwL(),-Wk*(Q*m)*row);      
      }
    } // end t has solid
    
    dQsf+=JacCol(v.ml(),-0.5*dTwdx*Q);
    dQsf+=JacCol(v.mr(),-0.5*dTwdx*Q);    
    return dQsf;
  }

  vc HopkinsHeatSource::HeatTransferCoefQ(const Cell& v) const{
    TRACE(10,"HopkinsHeatSource::HeatTransferCoefQ(const Cell& v)");
    const us& Nf=v.gc->Nf();

    // TRACE(100,"dTwdx:"<<dTwdx);
    const d& rh=v.vrh;    
    const d T0=v.T()(0);
    const d Pr0=v.gc->gas().pr(T0);
    const d cp0=v.gc->gas().cp(T0);    
    const d p0=v.p()(0)+v.gc->p0();
    const d rho0=v.gc->gas().rho(T0,p0);
    const d kappa0=v.gc->gas().kappa(T0);
    const d mu0=v.gc->gas().mu(T0);
    vc htcoefQ(Nf+1,fillwith::zeros);

    htcoefQ(0)=cp0*zeroheatQ;
    if(Nf>0){
      vd omgvec=v.gc->getomg()*linspace(1,Nf,Nf);
      vd deltak=sqrt(2*kappa0/(rho0*cp0*omgvec));
      vd deltanu=sqrt(2*mu0/(rho0*omgvec));
      vc fnu=rf.fx(rh/deltanu); // Viscous rott function
      vc fk=rf.fx(rh/deltak); // Thermal rott function
      htcoefQ.subvec(1,Nf)=cp0*(1/(1-Pr0))*(fnu/(1-fnu)-Pr0*fk/(1-fk));
    }
    // No time-average part here.
    return htcoefQ;
  }
  vc HopkinsHeatSource::HeatTransferCoefH(const Cell& v) const{
    TRACE(8,"HopkinsHeatSource::HeatTransferCoefH()");
    assert(zeroheatH_funptr);
    const us& Nf=v.gc->Nf();
    d T0=v.T()(0);	// Time-averaged temperature
    d kappa0=v.gc->gas().kappa(T0);
    d p0=v.p()(0)+v.gc->p0();
    d rho0=v.gc->gas().rho(T0,p0);
    d cp0=v.gc->gas().cp(T0);
    const d& rh=v.vrh;

    d omg=v.gc->getomg();
    vc htcoefH(Nf+1,fillwith::zeros);
    // Checked, this is correct
    htcoefH(0)=(*zeroheatH_funptr)(kappa0,rh);
    if(Nf>0){
      vd omgvec=omg*linspace(1,Nf,Nf);
      vd deltak=sqrt(2*kappa0/(rho0*cp0*omgvec));
      vc fk=rf.fx(rh/deltak); // Thermal rott function
      htcoefH.subvec(1,Nf)=I*rho0*cp0*omgvec%(fk/(1-fk));
    }
    // TRACE(100,"htcoefH0:"<<htcoefH(0));
    return htcoefH;
  }

  
} // namespace tube

//////////////////////////////////////////////////////////////////////
