// solidenergy.h
//
// Author: J.A. de Jong 
//
// Description:
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef SOLIDENERGY_H
#define SOLIDENERGY_H
#include "tubeequation.h"
#include "constants.h"

namespace tube{

  class SolidTPrescribed:public Equation
  {
    d Tsmirror=constants::T0;
  public:
    SolidTPrescribed(const Cell& v):Equation(v){}
    virtual void init();
    virtual enum EqType getType() const { return EqType::Sol;}
    tasystem::JacRow jac() const;
    
    vd error() const;			// Error in Solidenergy equation at
                                // node i
    vd extrapolateHeatFlow() const;
    tasystem::JacRow dExtrapolateHeatFlow() const;
    void setTs(d Ts) {Tsmirror=Ts;}
    vd kappaL() const;
    vd kappaR() const;
    void show() const;
  private:
    tasystem::JacCol dTsi() const;
    
  };
}

#endif // SOLIDENERGY_H
//////////////////////////////////////////////////////////////////////
