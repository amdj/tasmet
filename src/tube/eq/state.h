#pragma once
// File stateeq.h
#include "tubeequation.h"


namespace tube{
  SPOILNAMESPACE

  // class State:public TubeEquation
  // {
  // public:
  //   virtual void init(const WeightFactors&,const Tube&);
  //   tasystem::JacRow jac() const;
  //   vd error() const;			// Error in momentum equation at node i
  // private:
  //   tasystem::JacCol drhoi() const;
  //   virtual tasystem::JacCol dpL() const;
  //   virtual tasystem::JacCol dpR() const;    
  //   tasystem::JacCol dTi() const;
  // };

  class State:public TubeEquation
  {
    d WLi=0,WLim1=0; 
  public:
    State(const TubeVertex& v):TubeEquation(v){}
    virtual void init();   
    tasystem::JacRow jac() const;
    vd error() const;
  private:
    tasystem::JacCol dpL() const;
    tasystem::JacCol drhoi() const;
    tasystem::JacCol dTi() const;
    tasystem::JacCol drhoim1() const;
    tasystem::JacCol dTim1() const;
    tasystem::JacCol drhoip1() const;
    tasystem::JacCol dTip1() const;
  };
  class StateR:public TubeEquation
  {
  public:
    StateR(const TubeVertex& v):TubeEquation(v){}
    tasystem::JacRow jac() const;
    vd error() const;			// Error in 
                                // i
    virtual void init();   
  protected:
    tasystem::JacCol dpR() const;    
    tasystem::JacCol drhoR() const;
    tasystem::JacCol dTR() const;

  };
  
}