#pragma once
#include "tubeequation.h"


namespace tube{
  class Solidenergy:public Equation
  {
  public:
    Solidenergy(Tube* tube,TubeVertex* gp);
    ~Solidenergy();
    dmat operator()();
    vd Error();			// Error in Solidenergy equation at node i
    dmat dTsi();    
    
  };
}
