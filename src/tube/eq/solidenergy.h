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

namespace solids{
  class Solid;
}

namespace tube{
  
  class SolidEnergy:public Equation
  {
    const Tube* t=nullptr;
    d ksfrac;
    const solids::Solid* solid;
  public:
    SolidEnergy(const Cell& v,const solids::Solid* s,d ksfrac=1.0);
    virtual void init();
    virtual tasystem::JacRow jac() const;
    virtual enum EqType getType() const { return EqType::Sol;}
    virtual void show() const; 
    virtual vd error() const;			// Error in Energy equation at node i
    virtual void domg(vd&) const;

    static vd extrapolateHeatFlow(const Cell&); 
    static tasystem::JacRow dExtrapolateHeatFlow(const Cell&);
    vd kappaLt(const Cell&) const;
    vd kappaRt(const Cell&) const;

    vd QL(const Cell&) const;              // Heat conduction trough left wall
    vd QR(const Cell&) const;              // Heat conduction trough right wall
    tasystem::JacRow dQL(const Cell&) const;
    tasystem::JacRow dQR(const Cell&) const;

  };				// SolidEnergy

}

#endif // SOLIDENERGY_H
//////////////////////////////////////////////////////////////////////
