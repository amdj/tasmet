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
    d WgRr=0,WgRl=0,WgLr=0,WgLl=0;
    
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
  private:
    // tasystem::JacRow dQR() const;
    // tasystem::JacRow dQL() const;

    tasystem::JacCol dUi() const;
    tasystem::JacCol dTi() const;
    tasystem::JacCol drhoi() const;

    tasystem::JacCol dpR() const;
    tasystem::JacCol dpL() const;
    
    tasystem::JacCol drhoR() const;
    tasystem::JacCol dUR() const;
    tasystem::JacCol dTR() const;

    tasystem::JacCol drhoL() const;
    tasystem::JacCol dUL() const;
    tasystem::JacCol dTL() const;

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

    vd QR() const;              // Heat conduction trough right wall
    vd QL() const;              // Heat conduction trough right wall

  };
}      // namespace tube
#endif /* _ENERGYEQ_H_ */
