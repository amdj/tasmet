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
    Equation(const tasystem::Globalconf& gc);
    virtual ~Equation() {}

    virtual vd Error()=0;
    virtual dmat Jac()=0;		// Returns the local Jacobian of this equation
    const tasystem::Globalconf& gc;    
    const us& Ns;
  };

  
}




