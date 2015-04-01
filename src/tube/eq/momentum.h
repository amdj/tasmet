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

namespace tube{
  SPOILNAMESPACE
  class Cell;
  class Tube;
  class DragResistance;

  class Momentum:public Equation
  {
    const DragResistance* drag=nullptr;
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

  };
}

#endif // MOMENTUM_H
//////////////////////////////////////////////////////////////////////
