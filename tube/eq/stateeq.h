#pragma once
// File stateeq.h
#include "tubeequation.h"


namespace tube{
  SPOILNAMESPACE

  // class State:public TubeEquation
  // {
  // public:
  //   virtual void init(const WeightFactors&,const Tube&);
  //   tasystem::JacRow jac(const TubeVertex&) const;
  //   vd error(const TubeVertex&) const;			// Error in momentum equation at node i
  // private:
  //   tasystem::JacCol drhoi(const TubeVertex&) const;
  //   virtual tasystem::JacCol dpL(const TubeVertex&) const;
  //   virtual tasystem::JacCol dpR(const TubeVertex&) const;    
  //   tasystem::JacCol dTi(const TubeVertex&) const;
  // };

  class StateL:public TubeEquation
  {
    d WLi=0,WLip1=0,WLim1=0; 
  public:
    StateL(const TubeVertex& v):TubeEquation(v){}
    virtual void init(const WeightFactors&,const Tube&);   
    tasystem::JacRow jac(const TubeVertex&) const;
    vd error(const TubeVertex&) const;			// Error in momentum equation at node i
  private:
    tasystem::JacCol dpL(const TubeVertex&) const;
    tasystem::JacCol drhoi(const TubeVertex&) const;
    tasystem::JacCol dTi(const TubeVertex&) const;
    tasystem::JacCol drhoim1(const TubeVertex&) const;
    tasystem::JacCol dTim1(const TubeVertex&) const;
    tasystem::JacCol drhoip1(const TubeVertex&) const;
    tasystem::JacCol dTip1(const TubeVertex&) const;
  };
  class StateR:public TubeEquation
  {
  public:
    tasystem::JacRow jac(const TubeVertex&) const;
    vd error(const TubeVertex&) const;			// Error in momentum equation at node i
  protected:
    tasystem::JacCol drhoi(const TubeVertex&) const;
    virtual tasystem::JacCol dpR(const TubeVertex&) const;    
    tasystem::JacCol dTi(const TubeVertex&) const;
    tasystem::JacCol drhoim1(const TubeVertex&) const;
    tasystem::JacCol dTim1(const TubeVertex&) const;

  };
  
}
