#include "isentropiceq.h"
#include "tubevertex.h"



namespace tube{
  void Isentropic::init(const Tube& t) {
    TRACE(6,"Isentropic::init(t)");
  }
  JacRow Isentropic::jac(const TubeVertex& v) const{
    TRACE(6,"Isentropic::jac()");
    JacRow jac(v.p,2);
    jac+=dpi(v);
    jac+=drhoi(v);
    return jac;
  }
  vd Isentropic::error(const TubeVertex& v) const {
    TRACE(6,"Isentropic::Error()");
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
  JacCol Isentropic::dpi(const TubeVertex& v) const {
    TRACE(1,"Isentropic::dpi()");
    JacCol dpi(v.p,arma::eye(v.gc->Ns,v.gc->Ns));
    d p0=v.gc->p0;    
    dpi.data()*=1/p0;
    return dpi;
  }
  JacCol Isentropic::drhoi(const TubeVertex& v) const {
    TRACE(1,"Isentropic::drhoi()"); 
    JacCol drhoi(v.rho);

    d T0=v.gc->T0;
    d p0=v.gc->p0;
    d rho0=v.gc->gas.rho(T0,p0);
    d gamma=v.gc->gas.gamma(T0);
    drhoi+=-1.0*(gamma/rho0)*v.gc->fDFT*
      diagmat(pow(v.rho.tdata()/rho0,(gamma-1.0)))*v.gc->iDFT;
    return drhoi;
  }

} // namespace tube
