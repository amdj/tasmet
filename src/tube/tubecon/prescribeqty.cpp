#include "prescribeqty.h"
#include "jacrow.h"
#include "var.h"

namespace tube{

  using tasystem::JacRow;
  using tasystem::JacCol;
  using variable::var;

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

  void PrescribeddxQty::set(us eqnr,const var& Qi,const var& Qj,const var& Qk,\
                            d xi,d xj,d xk,const var& vals){
    TRACE(15,"PrescribeQty::set()");

    this->Qi=&Qi;
    this->Qj=&Qj;
    this->Qk=&Qk;

    d dxj=xj-xi;
    d dxk=xk-xi;

    d dxj_dxk=dxj/dxk;
    d dxj_dxk_sq=pow(dxj_dxk,2);

    d denom=dxj*(1-dxj_dxk);

    Wi=(dxj_dxk_sq-1)/denom;
    Wj=1/denom;
    Wk=-dxj_dxk_sq/denom;

    this->eqnr=eqnr;
    this->vals=vals;
  }

  vd PrescribeddxQty::error() const {
    TRACE(5,"PrescribeQty::error()");
    return Wi*Qi->adata()+Wk*Qk->adata()+Wj*Qj->adata()-vals();
  }
  JacRow PrescribeddxQty::jac() const {
    TRACE(5,"PrescribeQty::jac()");
    JacRow jac(eqnr,3);
    jac+=JacCol(*Qi,Wi*eye<dmat>(Qi->gc().Ns(),Qi->gc().Ns()));
    jac+=JacCol(*Qj,Wj*eye<dmat>(Qi->gc().Ns(),Qi->gc().Ns()));
    jac+=JacCol(*Qk,Wk*eye<dmat>(Qi->gc().Ns(),Qi->gc().Ns()));
    return jac;
  }
} // namespace tube







