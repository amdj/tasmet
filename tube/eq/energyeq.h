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
  public:
    d Htot(const TubeVertex&) const;             // Total enthalpy flow through this node
    virtual void init(const WeightFactors&);
    virtual JacRow jac(const TubeVertex&) const;
    virtual enum EqType getType() const { return EqType::Ene;}
    void show() const; 
    vd error(const TubeVertex&) const;			// Error in Energy equation at node i
    virtual void domg(const TubeVertex&,vd&) const;

    JacCol dUi(const TubeVertex&) const;
    JacCol dTi(const TubeVertex&) const;
    JacCol drhoi(const TubeVertex&) const;

    virtual JacCol dpR(const TubeVertex&) const;
    virtual JacCol dpL(const TubeVertex&) const;
    
    JacCol drhoip1(const TubeVertex&) const;
    JacCol dUip1(const TubeVertex&) const;
    JacCol dTip1(const TubeVertex&) const;

    JacCol drhoim1(const TubeVertex&) const;
    JacCol dUim1(const TubeVertex&) const;
    JacCol dTim1(const TubeVertex&) const;
    
    vd kappaL(const TubeVertex&) const;		// Thermal conducticity at the left
				// boundary of the vertex
    vd kappaR(const TubeVertex&) const;		// Thermal conductivity at the right
				// boundary of the vertex
    d gamma(const TubeVertex&) const;			// Time-avg ratio of specific heats
  };
}      // namespace tube
#endif /* _ENERGYEQ_H_ */
