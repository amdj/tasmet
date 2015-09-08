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
  class Cell;

  class  PrescribeQty{
    const tasystem::var* toprescribe;
    tasystem::var vals;
    us eqnr;
  public:
    PrescribeQty(){}
    // Vals: variables to coerce to
    PrescribeQty(const tasystem::var& vals){setVals(vals);} // 
    // toprescribe: a reference to the variable which has to be prescribed
    void set(us eqnr,const tasystem::var& toprescribe,const tasystem::var& vals);
    void set(us eqnr,const tasystem::var& toprescribe);
    void setVals(const tasystem::var& vals);
    vd error() const;
    tasystem::JacRow jac() const;
    void updateNf() {vals.updateNf();}
    void setGc(const tasystem::Globalconf& gc){vals.setGc(gc);}
    const tasystem::var& getVals() const {return vals;}
  };

} // namespace tube

#endif /* _PRESCRIBEQTY_H_ */
