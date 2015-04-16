// momentum.h
//
// Author: J.A. de Jong 
//
// Description:
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef MOMENTUM_H
#define MOMENTUM_H
#include "tubeequation.h"

// Define this variable to disable all drag
// #define NODRAG

namespace tube{
  SPOILNAMESPACE
  class Cell;
  class Tube;
  class DragResistance;

  class Momentum:public Equation
  {
    const Tube* t=nullptr;
  public:
    d Wddt=0,Wpi=0,Wpim1=0;

    Momentum(const Cell& v):Equation(v){}
    ~Momentum();
    virtual void init();
    virtual tasystem::JacRow jac() const;
    virtual enum EqType getType() const { return EqType::Mom;}    
    virtual vd error() const;			// Error in momentum equation at node i
    virtual void show() const;
    virtual void domg(vd& domg) const;

    static vd extrapolatePressure(const Cell&);
    static tasystem::JacRow dExtrapolatePressure(const Cell&);
    // Put in leftmost or rightmost cell and you obtain an
    // extrapolation of the momentumflow to the closest side
    static vd extrapolateMomentumFlow(const Cell&);
    static tasystem::JacRow dExtrapolateMomentumFlow(const Cell&);
  };
}

#endif // MOMENTUM_H
//////////////////////////////////////////////////////////////////////
