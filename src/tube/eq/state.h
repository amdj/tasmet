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
    d Wl=0,Wr=0; 
  public:
    State(const TubeVertex& v):TubeEquation(v){}
    virtual void init();   
    tasystem::JacRow jac() const;
    vd error() const;
  };
  class StateR:public TubeEquation
  {
  public:
    StateR(const TubeVertex& v):TubeEquation(v){}
    tasystem::JacRow jac() const;
    vd error() const;			// Error in 
                                // i
    virtual void init();   
  };
  
}
