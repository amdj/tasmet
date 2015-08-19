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

  class BcCell;
  
  class State:public Equation
  {
  public:
    State(const Cell& v):Equation(v){}
    virtual void init(){}
    void show() const;
    tasystem::JacRow jac() const;
    vd error() const;
    EqType getType() const {return EqType::Sta;}
  };
  class StateBc:public Equation
  {
    const BcCell& v;
  public:
    StateBc(const BcCell& v);
    virtual void init(){}
    void show() const;
    tasystem::JacRow jac() const;
    vd error() const;
    EqType getType() const {return EqType::BcEqStateBc;}
  };
  
}
#endif // STATE_H
//////////////////////////////////////////////////////////////////////
