#include "isentropiceq.h"
#include "tubevertex.h"



namespace tube{
  void Isentropic::init(const Tube& t) {
    TRACE(6,"Isentropic::init(t)");
  }
  JacRow Isentropic::jac(const TubeVertex& v) const{
    TRACE(6,"Isentropic::jac()");
    JacRow jac(dofnr,2);
    TRACE(0,"Isentropic, dofnr jac:"<< dofnr);
    jac+=dpL(v);
    jac+=dpR(v);    
    jac+=drhoi(v);
    return jac;
  }
  vd Isentropic::error(const TubeVertex& v) const {
    TRACE(6,"Isentropic::Error()");
    vd err(v.gc->Ns(),fillwith::zeros);
    d T0=v.gc->T0;
    d p0=v.gc->p0;
    d rho0=v.gc->gas.rho(T0,p0);
    d gamma=v.gc->gas.gamma(T0);
    vd p0vec_freqdomain=vd(v.gc->Ns(),fillwith::zeros);
    p0vec_freqdomain(0)=v.gc->p0;

    // Integrated form
    // err+=(v.eWisrho*p0vec_freqdomain+(v.eWispL*v.pL()()+v.eWispR*v.pR()()))/p0;
    // err+=-v.eWisrho*v.gc->fDFT*pow(v.rho.tdata()/rho0,gamma);

    err+=(p0vec_freqdomain+0.5*(v.pL()()+v.pR()()))/p0;
    err+=-v.gc->fDFT*pow(v.rho.tdata()/rho0,gamma);

    
    TRACE(6,"Isentropic::Error() done");
    return err;
  }
  JacCol Isentropic::dpL(const TubeVertex& v) const {
    TRACE(1,"Isentropic::dpi()");
    JacCol dpL(v.pL(),arma::eye(v.gc->Ns(),v.gc->Ns()));
    d p0=v.gc->p0;    
    // Integrated form
    // dpL.data()*=v.eWispL/p0;
    dpL.data()*=0.5/p0;
    return dpL;
  }
  JacCol Isentropic::dpR(const TubeVertex& v) const {
    TRACE(1,"Isentropic::dpi()");
    JacCol dpR(v.pR(),arma::eye(v.gc->Ns(),v.gc->Ns()));
    d p0=v.gc->p0;    
    // Integrated form
    // dpR.data()*=v.eWispR/p0;
    dpR.data()*=0.5/p0;
    return dpR;
  }  JacCol Isentropic::drhoi(const TubeVertex& v) const {
    TRACE(1,"Isentropic::drhoi()"); 
    JacCol drhoi(v.rho);

    d T0=v.gc->T0;
    d p0=v.gc->p0;
    d rho0=v.gc->gas.rho(T0,p0);
    d gamma=v.gc->gas.gamma(T0);
    // Integrated form
    // drhoi+=-1.0*v.eWisrho*(gamma/rho0)*v.gc->fDFT*
      // diagmat(pow(v.rho.tdata()/rho0,(gamma-1.0)))*v.gc->iDFT;
    drhoi+=-1.0*(gamma/rho0)*v.gc->fDFT*
      diagmat(pow(v.rho.tdata()/rho0,(gamma-1.0)))*v.gc->iDFT;
    return drhoi;
  }

} // namespace tube
