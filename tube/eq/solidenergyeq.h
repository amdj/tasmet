#pragma once
#include "tubeequation.h"


namespace tube{
  class Solidenergy:public TubeEquation
  {
  public:
    Solidenergy();
    ~Solidenergy();
    vd error(const TubeVertex&) const;			// Error in Solidenergy equation at node i
    dmat dTsi(const TubeVertex&) const;    
  };
}
