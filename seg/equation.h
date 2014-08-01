#pragma once
#include <vtypes.h>
#include "segbase.h"

// The only thing an equation does is return a Jacobian, and we are
// able to compute the error. That's it. And we only promise that
// these quantities can be obtained.

namespace segment{
  SPOILNAMESPACE

  class Equation
  {
  public:
    virtual void show() const { cout << "Empty equation description. From equation.h.\n";}
  };

  
}




