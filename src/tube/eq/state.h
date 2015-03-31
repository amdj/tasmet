// state.h
//
// Author: J.A. de Jong 
//
// Description:
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef STATE_H
#define STATE_H
#include "tubeequation.h"


namespace tube{
  SPOILNAMESPACE

  class State:public Equation
  {
  public:
    State(const Cell& v):Equation(v){}
    virtual void init(){}
    tasystem::JacRow jac() const;
    vd error() const;
    EqType getType() const {return EqType::Sta;}
  };
  
}
#endif // STATE_H
//////////////////////////////////////////////////////////////////////
