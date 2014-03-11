#pragma once
#include "../var/var.h"
#include <math_common.h>
#include "geom.h"

namespace tube{
  class Tube;
  class TubeVertex;

  class Equation{
  public:
    Equation(Tube* tube,TubeVertex* gp);
    ~Equation();
    d Sf,dxp,dxm;		// Fluid cs-area, dx+,dx-
    us i; 			// Current node
    Tube* tube;
    TubeVertex* vertex;		// Reference to parent (current gridpoint)
    const variable::varoperations& vop;
    const us& Ns;
    const Geom& geom;
    virtual vd Error()=0;
    virtual dmat operator()();		// Returns the local Jacobian of this equation
    dmat zero;			// Zeros matrix of right size
    virtual dmat drhoim1();	// Derivative of current equation to density at node i-1
    virtual dmat dUim1();	// Etc
    virtual dmat dTim1();
    virtual dmat dpim1();
    virtual dmat dTsim1();

    virtual dmat drhoi();	// Derivative of current equation to density at node i
    virtual dmat dUi();		// Etc
    virtual dmat dTi();
    virtual dmat dpi();
    virtual dmat dTsi();

    virtual dmat drhoip1();	// Derivative of current equation to pressure
    virtual dmat dUip1();
    virtual dmat dTip1();
    virtual dmat dpip1();
    virtual dmat dTsip1();
  };				// class Equation

} // namespace tube
