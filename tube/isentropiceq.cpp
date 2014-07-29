#include "isentropiceq.h"
#include "tubevertex.h"



namespace tube{
  Isentropic::Isentropic(TubeVertex& gp):TubeEquation(gp){
  }
  Isentropic::~Isentropic(){}
  vd Isentropic::Error(){
    TRACE(6,"Isentropic::Error()");
    TRACE(6,"Isentropic::Error() done");
    vd err(v.gc->Ns,fillwith::zeros);
    d T0=v.gc->T0;
    d p0=v.gc->p0;
    d rho0=v.gc->gas.rho(T0,p0);
    d gamma=v.gc->gas.gamma(T0);
    vd p0vec_freqdomain=vd(v.gc->Ns,fillwith::zeros);
    p0vec_freqdomain(0)=v.gc->p0;
    err+=(p0vec_freqdomain+v.p())/p0;
    err+=-1.0*v.gc->fDFT*pow(v.rho.tdata()/rho0,gamma);
    TRACE(6,"Isentropic::Error() done");
    return err;
  }
  dmat Isentropic::dpi(){
    TRACE(1,"Isentropic::dpi()");
    dmat dpi(v.gc->Ns,v.gc->Ns,fillwith::eye);
    d p0=v.gc->p0;    
    dpi=dpi/p0;
    return dpi;
  }
  dmat Isentropic::drhoi(){
    TRACE(1,"Isentropic::drhoi()"); 
    dmat drhoi(v.gc->Ns,v.gc->Ns,fillwith::zeros);

    d T0=v.gc->T0;
    d p0=v.gc->p0;
    d rho0=v.gc->gas.rho(T0,p0);

    d gamma=v.gc->gas.gamma(T0);
    drhoi+=-1.0*(gamma/rho0)*v.gc->fDFT*
      diagmat(pow(v.rho.tdata()/rho0,(gamma-1.0)))*v.gc->iDFT;
    return drhoi;
  }

} // namespace tube
