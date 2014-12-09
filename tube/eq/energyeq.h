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
    d WgiRr=0,WgiRl=0,WgiLr=0,WgiLl=0;
    
    d WcLl=0,WcLr=0,WcRl=0,WcRr=0; // Conduction weight factors

    d Wkin=0,WkinL=0,WkinR=0;
  public:
    Energy(const TubeVertex& v):TubeEquation(v){}
    d Htot() const;             // Total enthalpy flow through this node
    virtual void init();
    virtual tasystem::JacRow jac() const;
    virtual enum EqType getType() const { return EqType::Ene;}
    virtual void show() const; 
    virtual vd error() const;			// Error in Energy equation at node i
    virtual void domg(vd&) const;

    vd EkinR() const;
    vd EkinL() const;
    vd hR() const;              // Static enthalpy flowing trhough
                                // right cell wall
    vd hL() const;              // Static enthalpy flowing trhough
                                // left cell wall
    vd QR() const;              // Heat conduction trough right wall
    vd QL() const;              // Heat conduction trough right wall
    JacRow dQR() const;
    JacRow dQL() const;


    vd kappaL() const;		// Thermal conducticity at the left
    // boundary of the vertex
    vd kappaR() const;		// Thermal conductivity at the right
    // boundary of the vertex
    d gamma() const;			// Time-avg ratio of specific heats
  private:
    tasystem::JacCol dUi() const;
    tasystem::JacCol dTi() const;
    tasystem::JacCol drhoi() const;

    tasystem::JacCol dpR() const;
    tasystem::JacCol dpL() const;
    
    tasystem::JacCol drhoip1() const;
    tasystem::JacCol dUip1() const;
    tasystem::JacCol dTip1() const;

    tasystem::JacCol drhoim1() const;
    tasystem::JacCol dUim1() const;
    tasystem::JacCol dTim1() const;

  };
}      // namespace tube
#endif /* _ENERGYEQ_H_ */
