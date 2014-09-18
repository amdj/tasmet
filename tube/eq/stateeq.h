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
    virtual JacCol dpL(const TubeVertex&) const;
    virtual JacCol dpR(const TubeVertex&) const;    
    JacCol dTi(const TubeVertex&) const;
  };

  class StateL:public TubeEquation
  {
  public:
    virtual TubeEquation* copy() const {return new StateL(*this);}    
    JacRow jac(const TubeVertex&) const;
    vd error(const TubeVertex&) const;			// Error in momentum equation at node i
  private:
    JacCol dpL(const TubeVertex&) const;
    JacCol drhoi(const TubeVertex&) const;
    JacCol dTi(const TubeVertex&) const;
    JacCol drhoim1(const TubeVertex&) const;
    JacCol dTim1(const TubeVertex&) const;
    JacCol drhoip1(const TubeVertex&) const;
    JacCol dTip1(const TubeVertex&) const;
  };
  class StateR:public TubeEquation
  {
  public:
    virtual TubeEquation* copy() const {return new StateR(*this);}    
    JacRow jac(const TubeVertex&) const;
    vd error(const TubeVertex&) const;			// Error in momentum equation at node i
  private:
    JacCol drhoi(const TubeVertex&) const;
    JacCol dpR(const TubeVertex&) const;    
    JacCol dTi(const TubeVertex&) const;
    JacCol drhoim1(const TubeVertex&) const;
    JacCol dTim1(const TubeVertex&) const;

  };
  
}
