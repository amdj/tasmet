#pragma once
// File stateeq.h
#include "tubeequation.h"


namespace tube{
  class State:public Equation
  {
  public:
    State(Tube* tube,TubeVertex* gp);
    ~State();
    dmat operator()();
    vd Error();			// Error in momentum equation at node i
    dmat drhoi();
    dmat dpi();
    dmat dTi();
    
    
  };
}
