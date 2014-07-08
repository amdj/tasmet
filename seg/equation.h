#pragma once
#include "var.h"
#include <math_common.h>
#include "geom.h"

namespace segment{
  SPOILNAMESPACE
  class Vertex;

  class Equation
  {
  public:
    Equation(){}
    virtual ~Equation() {}

    virtual vd Error()=0;
    virtual dmat Jac()=0;		// Returns the local Jacobian of this equation
    const tasystem::Globalconf* gc=NULL;
    void Init(const tasystem::Globalconf& g){gc=&g;}
  };

  
}




