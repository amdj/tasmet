#include "cell.h"
#include "jacrow.h"
#include "isentropic.h"

#define iDFT (v.gc->iDFT())
#define fDFT (v.gc->fDFT())

namespace duct{
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;

  void Isentropic::init() {
    TRACE(6,"Isentropic::init(t)");
  }
  void Isentropic::show() const{
    cout << "----------------- Isentropic\n";
  }
 vd Isentropic::error() const {
    TRACE(6,"Isentropic::Error()");
    vd err(v.gc->Ns(),fillwith::zeros);
    d T0=v.gc->T0();
    d p0=v.gc->p0();
    d rho0=v.gc->gas().rho(T0,p0);
    d gamma=v.gc->gas().gamma(T0);

    err+=v.p()()/p0;
    err(0)+=1;
    err+=-fDFT*pow(v.rho().tdata()/rho0,gamma);
    TRACE(6,"Isentropic::Error() done");
    return err;
  }
  JacRow Isentropic::jac() const{
    TRACE(6,"Isentropic::jac()");
    JacRow jac(dofnr,3);
    TRACE(0,"Isentropic, dofnr jac:"<< dofnr);
    d p0=v.gc->p0();    
    jac+=JacCol(v.p(),eye()/p0);
    d T0=v.gc->T0();
    d rho0=v.gc->gas().rho(T0,p0);
    d gamma=v.gc->gas().gamma(T0);
    // Integrated form
    jac+=JacCol(v.rho(),-(gamma/rho0)*fDFT*
                 diagmat(pow(v.rho().tdata()/rho0,(gamma-1.0)))*iDFT);
    return jac;
  }

} // namespace duct


