#pragma once
// File stateeq.h
#include "tubeequation.h"


namespace tube{
  SPOILNAMESPACE
  class State:public TubeEquation
  {
  public:
    vd error(const TubeVertex&) const;			// Error in momentum equation at node i
    dmat drhoi(const TubeVertex&) const;
    dmat dpi(const TubeVertex&) const;
    dmat dTi(const TubeVertex&) const;
    
    
  };
}
