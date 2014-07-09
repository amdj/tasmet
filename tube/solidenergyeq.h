#pragma once
#include "tubeequation.h"


namespace tube{
  class Solidenergy:public TubeEquation
  {
  public:
    Solidenergy(TubeVertex& gp);
    ~Solidenergy();
    vd Error();			// Error in Solidenergy equation at node i
    dmat dTsi();    
    
  };
}
