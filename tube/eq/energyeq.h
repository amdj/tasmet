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

    // The following are for boundary conditions
    d WgUip1pL=0,WgUim1pR=0;
    
    d Wc1=0,Wc2=0,Wc3=0,Wc4=0; // Conduction weight factors
    d Wkini=0,Wkinim1=0,Wkinip1=0;
    d wLl=0,wRr=0,wLr=0,wRl=0;		// Basic weight functions
    d wL0=0,wL1=0,wRNm1=0,wRNm2=0;    	// Special boundary weight functions

  public:
    Energy(const TubeVertex& v):TubeEquation(v){}
    d Htot() const;             // Total enthalpy flow through this node
    virtual void init(const WeightFactors&,const Tube&);
    virtual tasystem::JacRow jac() const;
    virtual enum EqType getType() const { return EqType::Ene;}
    void show() const; 
    vd error() const;			// Error in Energy equation at node i
    virtual void domg(vd&) const;

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
    
    vd kappaL() const;		// Thermal conducticity at the left
				// boundary of the vertex
    vd kappaR() const;		// Thermal conductivity at the right
				// boundary of the vertex
    d gamma() const;			// Time-avg ratio of specific heats
  };
}      // namespace tube
#endif /* _ENERGYEQ_H_ */
