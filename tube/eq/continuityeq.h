#pragma once
#include "tubeequation.h"

namespace tube{



  class Continuity:public TubeEquation{	// Tube continuity equation 
  public:
    virtual vd error(const TubeVertex&) const;			// Error in this equation at this node
    virtual dmat drhoip1(const TubeVertex&) const; // Derivative of continuity equation to density at node
		    // i + 1
    virtual dmat drhoi(const TubeVertex&) const;	// Derivative of continuity equation to density at
			// novirtual de i
    virtual dmat dUi(const TubeVertex&) const;	// Derivative of continuity equation to Volume flow at i
		// node (const TubeVertex&only nonzero for nonconstant grids)
    virtual dmat drhoim1(const TubeVertex&) const;	// Derivative of continuity equation to density at
			// node i - 1
    virtual dmat dUip1(const TubeVertex&) const;	// Derivative of continuity equation to Volume flow
			// at node i + 1
    virtual dmat dUim1(const TubeVertex&) const;	// Derivative of continuity equation to Volume flow
			// at node i - 1
    virtual dmat drhoip2(const TubeVertex&) const;
    virtual dmat drhoim2(const TubeVertex&) const;

  };				// Continuity class
}				// Namespace tube

