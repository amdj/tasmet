// hopkinsheat.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
// 
//////////////////////////////////////////////////////////////////////


#include "hopkinsheat.h"
#include "cell.h"
#include "tube.h"
#include "geom.h"

/* For obtaining time-averaged heat transfer, we have to find the
   time-average of the following parameter:
 
   i * omg * f_x
   ------------
   (1 - f_x )


   Heat transfer is related to temperature difference, according to the
   Hopkins papers with:

   Qin = H*(Tw-T) + Q*dTwdx*u


*/



namespace tube{
  namespace H{
    SPOILNAMESPACE
    namespace{ 
      d zeroheat_vert(d kappa,d rh){
        return 3*kappa/pow(rh,2);
      }
      d zeroheat_circ(d kappa,d rh){
        return 2*kappa/pow(rh,2);
      }
      d zeroheat_inviscid(d,d){
        TRACE(2,"zeroheat_inviscid");
        return 0;
      }
    }
  } // namespace H

  
  
  using rottfuncs::RottFuncs;
  HopkinsHeatSource::HopkinsHeatSource(const Tube& t):
    cshape(t.geom().shape())
  {
    TRACE(10,"HopkinsHeatSource::HopKinsHeatSource(Tube)");
    setZeroFreq(cshape);
    if(t.geom().isBlApprox())
      rf=RottFuncs("blapprox");
    else
      rf=RottFuncs(t.geom().shape()); // Reinitialize thermoviscous functions with right shape
    dTwdx=zeros(t.geom().nCells());
    TRACE(11,"Exiting redefinition of Rottfuncs");
  }
  void HopkinsHeatSource::setZeroFreq(const string& shape){
    TRACE(10,"HopkinsHeatSource::setZeroFreq()");
    if(shape.compare("vert")==0){
      zeroheatH_funptr=&H::zeroheat_vert;
      zeroheatQ=0.2;      
    }
    else if(shape.compare("circ")==0){
      zeroheatH_funptr=&H::zeroheat_circ;
      zeroheatQ=1/3;      
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
  void HopkinsHeatSource::setdTwdx(const vd& dTwdx){
    this->dTwdx=dTwdx;
  }
  vd HopkinsHeatSource::heat(const Cell& v) const{
    TRACE(5,"HopkinsHeatSource::heat(v)");
    vd heat(v.gc->Ns(),fillwith::zeros);
    heat+=dTi(v)*(v.T()()-v.Ts()());
    heat+=dmi(v)*(0.5*(v.mL()()+v.mR()()) );    
    return heat;    
  }
  dmat HopkinsHeatSource::dTi(const Cell& v) const{
    TRACE(5,"HopkinsHeatSource::dTi(v)");
    variable::var htcoefH(*v.gc);
    htcoefH.setadata(HeatTransferCoefH(v));
    dmat dTi(v.gc->Ns(),v.gc->Ns(),fillwith::zeros);
    dTi=v.vVf*htcoefH.freqMultiplyMat();
    return dTi;
  }
  dmat HopkinsHeatSource::dmi(const Cell& v) const{
    TRACE(5,"HopkinsHeatSource::dUi(v)");
    variable::var htcoefQ(*v.gc);
    htcoefQ.setadata(HeatTransferCoefQ(v));
    dmat dmi(v.gc->Ns(),v.gc->Ns(),fillwith::zeros);
    dmi=-(v.xR-v.xL)*htcoefQ.freqMultiplyMat();
    return dmi;
  }  
  vc HopkinsHeatSource::HeatTransferCoefQ(const Cell& v) const{
    TRACE(10,"HopkinsHeatSource::HeatTransferCoefQ(const Cell& v)");
    const us& Nf=v.gc->Nf();

    // Obtain dTwdx
    d dTwdx=this->dTwdx(v.geti());
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

    htcoefQ(0)=cp0*zeroheatQ*dTwdx;
    if(Nf>0){
      vd omgvec=v.gc->getomg()*linspace(1,Nf,Nf);
      vd deltak=sqrt(2*kappa0/(rho0*cp0*omgvec));
      vd deltanu=sqrt(2*mu0/(rho0*omgvec));
      vc fnu=rf.fx(rh/deltanu); // Viscous rott function
      vc fk=rf.fx(rh/deltak); // Thermal rott function
      htcoefQ.subvec(1,Nf)=-cp0*(1/(1-Pr0))*(fnu/(1-fnu)-Pr0*fk/(1-fk))*dTwdx;
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
    htcoefH(0)=(*zeroheatH_funptr)(kappa0,rh)/rho0;
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
