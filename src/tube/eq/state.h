#pragma once
// File stateeq.h
#include "tubeequation.h"


namespace tube{
  SPOILNAMESPACE

  class State:public TubeEquation
  {
  public:
    State(const Cell& v):TubeEquation(v){}
    virtual void init(){}
    tasystem::JacRow jac() const;
    vd error() const;
  };
  
}
