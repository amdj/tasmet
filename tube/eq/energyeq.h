#pragma once
#ifndef _ENERGYEQ_H_
#define _ENERGYEQ_H_
#include "tubeequation.h"
#include "heat.h"

namespace tube{
  SPOILNAMESPACE
  class Energy:public TubeEquation
  {
    const HeatSource* heat;
  public:
    virtual TubeEquation* copy() const {return new Energy(*this);}    
    virtual void init(const Tube&);
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
