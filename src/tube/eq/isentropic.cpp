#include "cell.h"
#include "jacrow.h"
#include "isentropic.h"

#define iDFT (v.gc->iDFT)
#define fDFT (v.gc->fDFT)

namespace tube{
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;

  void Isentropic::init() {
    TRACE(6,"Isentropic::init(t)");
  }
  JacRow Isentropic::jac() const{
    TRACE(6,"Isentropic::jac()");
    JacRow jac(dofnr,3);
    TRACE(0,"Isentropic, dofnr jac:"<< dofnr);
    jac+=dp();    
    jac+=drho();
    return jac;
  }
  vd Isentropic::error() const {
    TRACE(6,"Isentropic::Error()");
    vd err(v.gc->Ns(),fillwith::zeros);
    d T0=v.gc->T0();
    d p0=v.gc->p0();
    d rho0=v.gc->gas().rho(T0,p0);
    d gamma=v.gc->gas().gamma(T0);
    vd p0vec_freqdomain=zeros();
    p0vec_freqdomain(0)=p0;

    // Integrated form
    // err+=(v.eWisrho*p0vec_freqdomain+(v.eWispL*v.pL()()+v.eWispR*v.pR()()))/p0;
    // err+=-v.eWisrho*fDFT*pow(v.rho.tdata()/rho0,gamma);

    err+=(p0vec_freqdomain+v.p()())/p0;
    err+=-fDFT*pow(v.rho().tdata()/rho0,gamma);
    
    TRACE(6,"Isentropic::Error() done");
    return err;
  }
  JacCol Isentropic::dp() const {
    TRACE(1,"Isentropic::dpi()");
    JacCol dp(v.p(),eye());
    d p0=v.gc->p0();    
    // Integrated form
    // dpL.data()*=v.eWispL/p0;
    dp.data()*=1/p0;
    return dp;
  }
  JacCol Isentropic::drho() const {
    TRACE(1,"Isentropic::drhoi()"); 
    JacCol drhoi(v.rho());
    d T0=v.gc->T0();
    d p0=v.gc->p0();
    d rho0=v.gc->gas().rho(T0,p0);
    d gamma=v.gc->gas().gamma(T0);
    // Integrated form
    // drhoi+=-1.0*v.eWisrho*(gamma/rho0)*fDFT*
      // diagmat(pow(v.rho.tdata()/rho0,(gamma-1.0)))*iDFT;
    drhoi+=-1.0*(gamma/rho0)*fDFT*
      diagmat(pow(v.rho().tdata()/rho0,(gamma-1.0)))*iDFT;
    return drhoi;
  }

} // namespace tube
