#pragma once
#ifndef _ENERGYEQ_H_
#define _ENERGYEQ_H_
#include "tubeequation.h"

namespace tube{
  SPOILNAMESPACE
  class HeatSource;
  
  class Energy:public Equation
  {
    const HeatSource* heat=nullptr;

  public:
    d Wddt=0,Wddtkin=0;
    d WRr=0,WRl=0,WLr=0,WLl=0;
    d WcLl=0,WcLr=0,WcRl=0,WcRr=0; // Conduction weight factors
    // d WkinLl=0,WkinLr=0,WkinRl=0,WkinRr=0;
  public:
    Energy(const Cell& v):Equation(v){}
    d Htot() const;             // Total enthalpy flow through this node
    virtual void init();
    virtual tasystem::JacRow jac() const;
    virtual enum EqType getType() const { return EqType::Ene;}
    virtual void show() const; 
    virtual vd error() const;			// Error in Energy equation at node i
    virtual void domg(vd&) const;

    vd kappaLt() const;		// Thermal conducticity at the left
    // boundary of the cell
    vd kappaRt() const;		// Thermal conductivity at the right
    // boundary of the cell
    d gamma() const;			// Time-avg ratio of specific heats

    vd extrapolateEnthalpyFlow() const;
    vd extrapolateHeatFlow() const;
    tasystem::JacRow dExtrapolateEnthalpyFlow() const;
    tasystem::JacRow dExtrapolateHeatFlow() const;
  private:

    // Total enthalpy flow through left cell wall
    vd mHL() const;

    // Total enthalpy flow through righ cell wall
    vd mHR() const;

    tasystem::JacRow dmHL() const;
    tasystem::JacRow dmHR() const;

    vd QL() const;              // Heat conduction trough left wall
    vd QR() const;              // Heat conduction trough right wall
    tasystem::JacRow dQL() const;
    tasystem::JacRow dQR() const;

    // tasystem::JacRow dQL() const;
    // tasystem::JacRow dQR() const;

  };
}      // namespace tube
#endif /* _ENERGYEQ_H_ */
