#pragma once
// File stateeq.h
#include "tubeequation.h"


namespace tube{
  SPOILNAMESPACE
  class State:public TubeEquation
  {
    JacRow jac(const TubeVertex&) const;
    vd error(const TubeVertex&) const;			// Error in momentum equation at node i
  private:
    JacCol drhoi(const TubeVertex&) const;
    JacCol dpi(const TubeVertex&) const;
    JacCol dTi(const TubeVertex&) const;
    
    
  };
}
