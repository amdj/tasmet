#pragma once
#ifndef _PRESCRIBEQTY_H_
#define _PRESCRIBEQTY_H_
#include "vtypes.h"
#include "var.h"

namespace tasystem{
  class JacRow;

}
namespace tube{
  SPOILNAMESPACE;
  class TubeVertex;

  class PrescribeQty{
    const variable::var* toprescribe;
    variable::var vals;
    us eqnr;
  public:
    PrescribeQty(){}
    PrescribeQty(const variable::var& vals){setVals(vals);}
    void set(us eqnr,const variable::var& toprescribe,const variable::var& vals);
    void set(us eqnr,const variable::var& toprescribe);
    void setVals(const variable::var& vals);
    vd error() const;
    tasystem::JacRow jac() const;
    void updateNf() {vals.updateNf();}
    void setGc(const tasystem::Globalconf& gc){vals.setGc(gc);}
    const variable::var& getVals() const {return vals;}
  };
  class PrescribeddxQty{
    const variable::var *Qi,*Qj,*Qk; // Quantities as position i,j,k
    d Wi,Wj,Wk;
    variable::var vals;
    us eqnr;
  public:
    PrescribeddxQty(){}
    #define varref const variable::var&
    void set(us eqnr,varref Qi,varref Qj,varref Qk,d xi,d xj, d xk,const variable::var& vals);
    #undef varref
    vd error() const;
    tasystem::JacRow jac() const;
    void updateNf() {vals.updateNf();}
    void setGc(const tasystem::Globalconf& gc){vals.setGc(gc);}
    const variable::var& getVals() const {return vals;}
  };

} // namespace tube

#endif /* _PRESCRIBEQTY_H_ */
