#include "hopkinsheat.h"
#include "tubevertex.h"



/* For obtaining time-averaged heat transfer, we have to find the
   time-average of the following parameter:
 
   i * omg * f_x
   ------------
   (1 - f_x )


   Heat transfer is related to temperature difference, according to the
   Hopkins papers with:

   Qin = H*(Tw-T) + Q*dTwdx*u


*/

namespace H{
  SPOILNAMESPACE
  d zeroheat_vert(d kappa,d rh){
    return 3*kappa/pow(rh,2);
  }
  d zeroheat_circ(d kappa,d rh){
    return 2*kappa/pow(rh,2);
  }
  d zeroheat_blapprox(d kappa,d rh){
    TRACE(2,"zeroheat_blapprox");
    return 2*kappa/pow(rh,2);
    // return 0;
  }
  d zeroheat_inviscid(d dummy,d dummy2){
    TRACE(2,"zeroheat_inviscid");
    return 0;
  }
}


namespace tube{
  using rottfuncs::RottFuncs;
  HopkinsHeatSource::HopkinsHeatSource(const Tube& t):
    cshape(t.geom.shape),
    rf(cshape)
  {
    TRACE(10,"HopkinsHeatSource::HopKinsHeatSource()");
    setZeroFreq(cshape);
  }
  HopkinsHeatSource::HopkinsHeatSource(const HopkinsHeatSource& o){
    cshape=o.cshape;
    rf=RottFuncs(cshape);
    setZeroFreq(cshape);
  }
  HopkinsHeatSource& HopkinsHeatSource::operator=(const HopkinsHeatSource& o){
    cshape=o.cshape;
    rf=RottFuncs(cshape);
    setZeroFreq(cshape);
    return *this;
  }
  void HopkinsHeatSource::setZeroFreq(const string& shape){
    TRACE(10,"HopkinsHeatSource::setZeroFreq()");
    if(shape.compare("vert")==0)
      zeroheatH_funptr=&H::zeroheat_vert;
    else if(shape.compare("circ")==0)
      zeroheatH_funptr=&H::zeroheat_circ;
    else if(shape.compare("blapprox")==0)
      zeroheatH_funptr=&H::zeroheat_blapprox;
    else if(shape.compare("inviscid")==0)
      zeroheatH_funptr=&H::zeroheat_inviscid;
    else
      {
        WARN("Warning: tube.geom.shape unknown for ZeroHeatH. Aborting...");
        abort();
      }
  }
  void HopkinsHeatSource::setdTwdx(const Geom& g,const vd& dTwdx){
    this->dTwdx=&dTwdx;
  }
  vd HopkinsHeatSource::heat(const TubeVertex& v) const{
    TRACE(5,"HopkinsHeatSource::heat(v)");
    vd heat(v.gc->Ns,fillwith::zeros);
    variable::var htcoefH(v.gc);
    variable::var htcoefQ(v.gc);
    htcoefH.set(HeatTransferCoefH(v));
    htcoefQ.set(HeatTransferCoefQ(v));    
    // TRACE(100,"TminTs:\n"<<v.T()-v.Ts());
    // if(v.i==0)
      // TRACE(25,"H freqmultiplymat:"<< htcoefQ.freqMultiplyMat());
    heat+=htcoefH.freqMultiplyMat()*(v.T()-v.Ts());
    heat+=htcoefQ.freqMultiplyMat()*(v.U()/v.lg.vSf);    
    return heat;    
  }
  dmat HopkinsHeatSource::dTi(const TubeVertex& v) const{
    TRACE(5,"HopkinsHeatSource::dTi(v)");
    variable::var htcoefH(v.gc);
    htcoefH.set(HeatTransferCoefH(v));
    dmat dTi(v.gc->Ns,v.gc->Ns,fillwith::zeros);
    dTi=htcoefH.freqMultiplyMat();
    return dTi;
  }
  dmat HopkinsHeatSource::dUi(const TubeVertex& v) const{
    TRACE(5,"HopkinsHeatSource::dUi(v)");
    variable::var htcoefQ(v.gc);
    htcoefQ.set(HeatTransferCoefQ(v));
    dmat dUi(v.gc->Ns,v.gc->Ns,fillwith::zeros);
    dUi=htcoefQ.freqMultiplyMat()/v.lg.vSf;
    return dUi;
  }  
  vc HopkinsHeatSource::HeatTransferCoefQ(const TubeVertex& v) const{
    const us& Nf=v.gc->Nf;
    vc htcoefQ(Nf+1,fillwith::zeros);

    // Obtain dTwdx
    d dTwdx=(*(this->dTwdx))(v.i);
    // TRACE(100,"dTwdx:"<<dTwdx);
    const d& rh=v.lg.vrh;    
    d T0=v.T(0);
    d Pr0=v.gc->gas.pr(T0);
    d cp0=v.gc->gas.cp(T0);    
    d p0=v.p(0)+v.gc->p0;
    d rho0=v.gc->gas.rho(T0,p0);
    d kappa0=v.gc->gas.kappa(T0);
    d mu0=v.gc->gas.mu(T0);
    if(Nf>0){
      vd omgvec=v.gc->getomg()*linspace(1,Nf,Nf);
      vd deltak=sqrt(2*kappa0/(rho0*cp0*omgvec));
      vd deltanu=sqrt(2*mu0/(rho0*omgvec));
      vc fnu=rf.fx(rh/deltanu); // Viscous rott function
      vc fk=rf.fx(rh/deltak); // Thermal rott function
      htcoefQ.subvec(1,Nf)=rho0*cp0*(1/(1-Pr0))*(fnu/(1-fnu)-Pr0*fk/(1-fk))*dTwdx;
    }
    // No time-average part here.
    return htcoefQ;
  }
  vc HopkinsHeatSource::HeatTransferCoefH(const TubeVertex& v) const{
    TRACE(8,"HopkinsHeatSource::HeatTransferCoefH()");
    const us& Nf=v.gc->Nf;
    d T0=v.T(0);	// Time-averaged temperature
    d kappa0=v.gc->gas.kappa(T0);
    d p0=v.p(0)+v.gc->p0;
    d rho0=v.gc->gas.rho(T0,p0);
    d cp0=v.gc->gas.cp(T0);
    const d& rh=v.lg.vrh;

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

  
}
