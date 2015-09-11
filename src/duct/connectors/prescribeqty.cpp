#include "prescribeqty.h"
#include "jacrow.h"
#include "var.h"

namespace duct{

  using tasystem::JacRow;
  using tasystem::JacCol;
  using tasystem::var;

  void PrescribeQty::set(us eqnr,const var& toprescribe,const var& vals){
    TRACE(15,"PrescribeQty::set()");
    set(eqnr,toprescribe);
    setVals(vals);
  }
  void PrescribeQty::set(us eqnr,const var& toprescribe){
    TRACE(15,"PrescribeQty::set()");
    this->toprescribe=&toprescribe;
    this->eqnr=eqnr;
    vals.setGc(toprescribe);
  }
  void PrescribeQty::setVals(const var& vals){
    this->vals=vals;
  }

  vd PrescribeQty::error() const {
    TRACE(5,"PrescribeQty::error()");
    return toprescribe->adata()-vals();
  }
  JacRow PrescribeQty::jac() const {
    TRACE(5,"PrescribeQty::jac()");
    JacRow jac(eqnr,1);
    jac+=JacCol(*toprescribe,eye<dmat>(toprescribe->gc().Ns(),toprescribe->gc().Ns()));
    return jac;
  }

} // namespace duct







