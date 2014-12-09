#include "tubevertex.h"
#include "hopkinslaminarduct.h"
#include "jacobian.h"
#include "solidenergy.h"

namespace tube{
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;


  void SolidTPrescribed::init() {
    TRACE(6,"SolidTPrescribed::init(t)");
    const Tube& t=v.getTube();
    if(t.getName().compare("HopkinsLaminarDuct")==0){
      const HopkinsLaminarDuct& d=dynamic_cast<const HopkinsLaminarDuct&>(t);
      Tsmirror=&d.Tmirror;
    }

    // Nope, we do nothing with weight functions
    
  }
  JacRow SolidTPrescribed::jac() const{
    TRACE(6,"SolidTPrescribed::jac()");
    JacRow jac(dofnr,1);
    TRACE(0,"Dofnr jac:"<< dofnr);
    jac+=dTsi();
    return jac;
  }
  vd SolidTPrescribed::error() const {		// Error in momentum equation
    TRACE(6,"SolidTPrescribed::Error()");
    vd error(v.gc->Ns(),fillwith::zeros);
    assert(v.gc!=NULL);
    error=v.Ts()();
    if(Tsmirror){
      error(0)-=(*Tsmirror)(v.geti());
    }
    else
      error(0)-=v.gc->T0;
    return error;
  }
  JacCol SolidTPrescribed::dTsi() const {
    TRACE(0,"SolidTPrescribed:dTsi()");
    // Set solid temperature to zero
    return JacCol(v.Ts(),arma::eye(v.gc->Ns(),v.gc->Ns()));
  }

} // namespace tube


