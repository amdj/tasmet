#pragma once
#include "tubeequation.h"

namespace tube{


  class Continuity:public Equation{	// Tube continuity equation 
  public:
    Continuity(const Tube& tube,const TubeVertex& gp);
    ~Continuity();
    //    dmat operator()();		// Link to Equation operator()()
    vd Error();			// Error in this equation at this node
    dmat drhoip1();	// Derivative of continuity equation to density at node i + 1
    dmat drhoi();	// Derivative of continuity equation to density at node i
    dmat dUi();	// Derivative of continuity equation to Volume flow at node i (only nonzero for nonconstant grids)
    dmat drhoim1();	// Derivative of continuity equation to density at node i - 1
    dmat dUip1();	// Derivative of continuity equation to Volume flow at node i + 1
    dmat dUim1();	// Derivative of continuity equation to Volume flow at node i - 1

  };				// Continuity class
}				// Namespace tube
