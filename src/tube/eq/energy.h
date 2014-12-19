#pragma once
#ifndef _ENERGYEQ_H_
#define _ENERGYEQ_H_
#include "tubeequation.h"
#include "heat.h"

namespace tube{
  SPOILNAMESPACE

  
  class Energy:public TubeEquation
  {
    const HeatSource* heat=NULL;

  public:
    d Wddt=0,Wddtkin=0;
    d WRr=0,WRl=0,WLr=0,WLl=0;
    d WcLl=0,WcLr=0,WcRl=0,WcRr=0; // Conduction weight factors
    d WkinLl=0,WkinLr=0,WkinRl=0,WkinRr;
  public:
    Energy(const TubeVertex& v):TubeEquation(v){}
    d Htot() const;             // Total enthalpy flow through this node
    virtual void init();
    virtual tasystem::JacRow jac() const;
    virtual enum EqType getType() const { return EqType::Ene;}
    virtual void show() const; 
    virtual vd error() const;			// Error in Energy equation at node i
    virtual void domg(vd&) const;

    vd kappaL() const;		// Thermal conducticity at the left
    // boundary of the vertex
    vd kappaR() const;		// Thermal conductivity at the right
    // boundary of the vertex
    d gamma() const;			// Time-avg ratio of specific heats

    vd extrapolateEnergyFlow() const;
    vd extrapolateHeatFlow() const;
    tasystem::JacRow dExtrapolateEnergyFlow() const;
    tasystem::JacRow dExtrapolateHeatFlow() const;
  private:
    vd ddtEtherm() const;
    tasystem::JacRow dddtEtherm() const;
    vd ddtEtot() const;
    tasystem::JacRow dddtEtot() const;

    vd EkinR() const;
    vd EkinL() const;

    vd hL() const;              // Static enthalpy flowing through
    // left cell wall
    vd hR() const;              // Static enthalpy flowing through
    // right cell wall
    vd HL() const;              // Total enthalpy flow through left
                                // cell wall
    vd HR() const;              // Total enthalpy flow through righ
                                // cell wall
    tasystem::JacRow dhL() const;
    tasystem::JacRow dhR() const;
    tasystem::JacRow dHL() const;
    tasystem::JacRow dHR() const;

    vd QL() const;              // Heat conduction trough left wall
    vd QR() const;              // Heat conduction trough right wall
    // tasystem::JacRow dQL() const;
    // tasystem::JacRow dQR() const;

  };
}      // namespace tube
#endif /* _ENERGYEQ_H_ */
