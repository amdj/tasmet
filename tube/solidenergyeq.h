#pragma once
#include "tubeequation.h"


namespace tube{
  class Solidenergy:public Equation
  {
  public:
    Solidenergy(const Tube& tube,const TubeVertex& gp);
    ~Solidenergy();
    dmat operator()();
    vd Error();			// Error in Solidenergy equation at node i
    dmat dTsi();    
    
  };
}
