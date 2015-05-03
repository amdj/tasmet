#pragma once
#ifndef _ENERGYEQ_H_
#define _ENERGYEQ_H_
#include "tubeequation.h"

namespace tube{
  SPOILNAMESPACE
  class HeatSource;
  
  class Energy:public Equation
  {
    const Tube* t=nullptr;
    d Wddt=0,Wddtkin=0,WheatQ;

  public:
    Energy(const Cell& v):Equation(v){}
    d Htot() const;             // Total enthalpy flow through this node
    virtual void init();
    virtual tasystem::JacRow jac() const;
    virtual enum EqType getType() const { return EqType::Ene;}
    virtual void show() const; 
    virtual vd error() const;			// Error in Energy equation at node i
    virtual void domg(vd&) const;

    static vd kappaLt(const Cell&);		// Thermal conducticity at the left
    // boundary of the cell
    static vd kappaRt(const Cell&);		// Thermal conductivity at the right
    // boundary of the cell
    d gamma() const;			// Time-avg ratio of specific heats

    static vd extrapolateEnthalpyFlow(const Cell&);
    static vd extrapolateHeatFlow(const Cell&); 
    static tasystem::JacRow dExtrapolateEnthalpyFlow(const Cell&);
    static tasystem::JacRow dExtrapolateHeatFlow(const Cell&);

    // Total enthalpy flow through left cell wall
    static vd mHL(const Cell&);
    // Total enthalpy flow through righ cell wall
    static vd mHR(const Cell&);
    static tasystem::JacRow dmHL(const Cell&);
    static tasystem::JacRow dmHR(const Cell&);

    static vd QL(const Cell&);              // Heat conduction trough left wall
    static vd QR(const Cell&);              // Heat conduction trough right wall
    static tasystem::JacRow dQL(const Cell&);
    static tasystem::JacRow dQR(const Cell&);
  private:


    // tasystem::JacRow dQL() const;
    // tasystem::JacRow dQR() const;

  };
}      // namespace tube
#endif /* _ENERGYEQ_H_ */
