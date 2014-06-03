#pragma once
#include "tubeequation.h"

namespace tube{


  class Continuity:public TubeEquation{	// Tube continuity equation 
  public:

    Continuity(const Tube& tube,TubeVertex& gp);
    ~Continuity();
  
    vd Error();			// Error in this equation at this node
    dmat drhoip1(); // Derivative of continuity equation to density at node
		    // i + 1
    dmat drhoi();	// Derivative of continuity equation to density at
			// node i
    dmat dUi();	// Derivative of continuity equation to Volume flow at i
		// node (only nonzero for nonconstant grids)
    dmat drhoim1();	// Derivative of continuity equation to density at
			// node i - 1
    dmat dUip1();	// Derivative of continuity equation to Volume flow
			// at node i + 1
    dmat dUim1();	// Derivative of continuity equation to Volume flow
			// at node i - 1
    dmat drhoip2();
    dmat drhoim2();

    d Wim1,Wi,Wip1;		// Weight functions for relative
				// contribution to the equations.
    d Wddt;  			// Weight function time-derivative term
  };				// Continuity class
}				// Namespace tube

