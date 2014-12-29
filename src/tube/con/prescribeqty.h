#pragma once
#ifndef _PRESCRIBEQTY_H_
#define _PRESCRIBEQTY_H_
#include "vtypes.h"

namespace variable{
  class var;
}
namespace tasystem{
  class JacRow;

}
namespace tube{
  SPOILNAMESPACE;
  class TubeVertex;

  class PrescribeQty{
    const variable::var* toprescribe;
    vd vals;
    us eqnr;
  public:
    PrescribeQty(){}
    void set(us eqnr,const variable::var& toprescribe,const vd& vals);
    vd error() const;
    tasystem::JacRow jac() const;
  };
  class PrescribeddxQty{
    const variable::var *Qi,*Qj,*Qk; // Quantities as position i,j,k
    d Wi,Wj,Wk;
    vd vals;
    us eqnr;
  public:
    PrescribeddxQty(){}
    #define varref const variable::var&
    void set(us eqnr,varref Qi,varref Qj,varref Qk,d xi,d xj, d xk,const vd& vals);
    #undef varref
    vd error() const;
    tasystem::JacRow jac() const;
  };

} // namespace tube

#endif /* _PRESCRIBEQTY_H_ */
