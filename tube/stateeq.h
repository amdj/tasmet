#pragma once
// File stateeq.h
#include "tubeequation.h"


namespace tube{
  SPOILNAMESPACE
  class State:public TubeEquation
  {
  public:
    State(TubeVertex& gp);
    ~State();

    vd Error();			// Error in momentum equation at node i
    dmat drhoi();
    dmat dpi();
    dmat dTi();
    
    
  };
}
