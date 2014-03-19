#pragma once
#include "../var/var.h"
#include <math_common.h>
#include "geom.h"

namespace tube{
  class Tube;
  class TubeVertex;

  class Equation{
  public:
    Equation(const Tube& tube,const TubeVertex& gp);
    ~Equation();
    Equation(const Equation& other); // Copy constructor
    
    dmat diagtmat(const variable::var& v); // Create diagonal matrix with time domain data from variable
    const us& i; 			// Current node
    const Tube& tube;
    const TubeVertex& vertex;		// Reference to parent (current gridpoint)
    const variable::varoperations& vop;
    const dmat& fDFT,iDFT,DDTfd;	// forward, backward dicrete fourier transform, derivative to time matrix (freq domain)
    const us& Ns;
    const Geom& geom;

    const d& vSf;			// Vertex fluid cross-sectional area
    const d& vSs;			// Vertex solid cross-sectional area
    const d& vVf;			// Vertex cell fluid volume
    const d& vVs;			// Vertex cell solid volume

    d SfR;			// Cross-sectional area of right face
    d SfL;			// Cross-sectional area of left  face
    const d& wLl;	// Weight functions (see docs for info)
    const d& wRr;	// Weight functions (see docs for info)
    const d& wLr;	// Weight functions (see docs for info)
    const d& wRl;	// Weight functions (see docs for info)
    
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
