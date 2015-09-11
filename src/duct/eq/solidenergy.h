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
#include "ductequation.h"
#include "constants.h"

namespace solids{
  class Solid;
}

namespace duct{
  class Duct;

  class SolidEnergy:public Equation
  {
    const Duct* t=nullptr;
    d ksfrac;			// Factor reducing
    const solids::Solid* solid;
    d Qin=0;
  public:
    SolidEnergy(const Cell& v,const solids::Solid* s,d ksfrac=1.0);
    virtual void init();
    virtual tasystem::JacRow jac() const;
    virtual enum EqType getType() const { return EqType::Sol;}
    void setQin(d Qin1) {Qin=Qin1;} 		// Time-averaged heat input in this
				// cell in Watts
    virtual void show() const; 
    virtual vd error() const;			// Error in Energy equation at node i
    virtual void domg(vd&) const;

    vd extrapolateHeatFlow() const; 
    tasystem::JacRow dExtrapolateHeatFlow() const;
    vd kappaLt() const;
    vd kappaRt() const;

    vd QL() const;              // Heat conduction trough left wall
    vd QR() const;              // Heat conduction trough right wall
    tasystem::JacRow dQL() const;
    tasystem::JacRow dQR() const;

  };				// SolidEnergy

}

#endif // SOLIDENERGY_H
//////////////////////////////////////////////////////////////////////
