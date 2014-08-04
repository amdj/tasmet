#pragma once
#ifndef _ENERGYEQ_H_
#define _ENERGYEQ_H_
#include "tubeequation.h"


namespace tube{
  SPOILNAMESPACE
  class Energy:public TubeEquation
  {
  public:
    virtual enum EqType getType() const { return EqType::Ene;}
    // void show() const; 
    vd error(const TubeVertex&) const;			// Error in Energy equation at node i
    virtual dmat dpi(const TubeVertex&) const;
    virtual dmat dUi(const TubeVertex&) const;
    virtual dmat dTi(const TubeVertex&) const;
    
    virtual dmat dpip1(const TubeVertex&) const;
    virtual dmat dpip2(const TubeVertex&) const;
    virtual dmat dpim2(const TubeVertex&) const;
    virtual dmat dUip1(const TubeVertex&) const;
    virtual dmat dTip1(const TubeVertex&) const;
    
    virtual dmat dUim1(const TubeVertex&) const;
    virtual dmat dpim1(const TubeVertex&) const;
    virtual dmat dTim1(const TubeVertex&) const;
    
    vd kappaL(const TubeVertex&) const;		// Thermal conducticity at the left
				// boundary of the vertex
    vd kappaR(const TubeVertex&) const;		// Thermal conductivity at the right
				// boundary of the vertex
    d gamma(const TubeVertex&) const;			// Time-avg ratio of specific heats
  };

}      // namespace tube
#endif /* _ENERGYEQ_H_ */
