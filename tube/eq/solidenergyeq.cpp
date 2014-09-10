#include "solidenergyeq.h"
#include "tubevertex.h"

namespace tube{

  SolidTPrescribed::SolidTPrescribed(){
    // Tset is false

  }
  void SolidTPrescribed::init(const Tube& t) {
    TRACE(6,"SolidTPrescribed::init(t)");

  }
  JacRow SolidTPrescribed::jac(const TubeVertex& v) const{
    TRACE(6,"SolidTPrescribed::jac()");
    JacRow jac(v.Ts,1);
    jac+=dTsi(v);
    return jac;
  }
  vd SolidTPrescribed::error(const TubeVertex& v) const {		// Error in momentum equation
    TRACE(6,"SolidTPrescribed::Error()");
    vd error(v.gc->Ns,fillwith::zeros);
    assert(v.gc!=NULL);
    error=v.Ts();
    if(Tsmirror.size()>0)
      error(0)-=Tsmirror(v.i);
    return error;
  }
  JacCol SolidTPrescribed::dTsi(const TubeVertex& v) const {
    TRACE(0,"SolidTPrescribed:dTsi()");
    // Set solid temperature to zero
    return JacCol(v.Ts,arma::eye(v.gc->Ns,v.gc->Ns));
  }

} // namespace tube


