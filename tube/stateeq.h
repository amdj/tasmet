#pragma once
// File stateeq.h
#include "tubeequation.h"


namespace tube{
  class State:public Equation
  {
  public:
    State(const Tube& tube,TubeVertex& gp);
    ~State();

    vd Error();			// Error in momentum equation at node i
    dmat drhoi();
    dmat dpi();
    dmat dTi();
    
    
  };
}
