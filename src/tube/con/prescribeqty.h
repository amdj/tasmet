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

} // namespace tube

#endif /* _PRESCRIBEQTY_H_ */
