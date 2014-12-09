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
    d Wgip1=0,Wgip=0,Wgim=0,Wgim1=0;

    
    d WcL=0,WcR=0,WcRR=0,WcLL=0; // Conduction weight factors

    d Wkini=0,Wkinim1=0,Wkinip1=0;
  public:
    Energy(const TubeVertex& v):TubeEquation(v){}
    d Htot() const;             // Total enthalpy flow through this node
    virtual void init(const WeightFactors&,const Tube&);
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
