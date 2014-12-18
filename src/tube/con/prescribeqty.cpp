#include "prescribeqty.h"
#include "jacobian.h"
#include "var.h"

namespace tube{

  using tasystem::JacRow;
  using tasystem::JacCol;
  using variable::var;

  void PrescribeQty::set(us eqnr,const var& toprescribe,const vd& vals){
    TRACE(15,"PrescribeQty::set()");

    this->toprescribe=&toprescribe;
    this->eqnr=eqnr;
    this->vals=vals;
  }

  vd PrescribeQty::error() const {
    TRACE(5,"PrescribeQty::error()");
    return toprescribe->adata()-vals;
  }
  JacRow PrescribeQty::jac() const {
    TRACE(5,"PrescribeQty::jac()");
    JacRow jac(eqnr,1);
    jac+=JacCol(*toprescribe,eye<dmat>(toprescribe->gc().Ns(),toprescribe->gc().Ns()));
    return jac;
  }
} // namespace tube
