#pragma once
// File stateeq.h
#include "tubeequation.h"


namespace tube{
  SPOILNAMESPACE
  class State:public TubeEquation
  {
  public:
    virtual TubeEquation* copy() const {return new State(*this);}    
    JacRow jac(const TubeVertex&) const;
    vd error(const TubeVertex&) const;			// Error in momentum equation at node i
  private:
    JacCol drhoi(const TubeVertex&) const;
    JacCol dpL(const TubeVertex&) const;
    JacCol dpR(const TubeVertex&) const;    
    JacCol dTi(const TubeVertex&) const;
    
    
  };
}
