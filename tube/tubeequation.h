#pragma once
#include "var.h"
#include <math_common.h>
#include "geom.h"

namespace segment{
  class Vertex;
}

namespace tube{
  SPOILNAMESPACE  

  class Tube;
  class TubeVertex;

  class Equation{
  public:
    Equation(const Tube& tube,TubeVertex& gp);
    ~Equation();
    Equation(const Equation& other); // Copy constructor
    dmat diagtmat(variable::var& v); // Create diagonal matrix with time domain data from variable
    const us& i; 			// Current node
    const Tube& tube;
    TubeVertex& vertex;		// Reference to parent (current gridpoint)

    segment::Vertex*& left;
    segment::Vertex*& right;
    
    const tasystem::Globalconf& gc;

    const dmat& fDFT,iDFT,DDTfd;	// forward, backward dicrete fourier transform, derivative to time matrix (freq domain)
    const us& Ns;
    
    vd getp0(); 		// Create a vector of zero-pressure data
    vd getp0t();   		// Same, but then time domain data
    const Geom& geom;
    const us& Ncells;		// Number of cells
    d vSf;			// Vertex fluid cross-sectional area
    d vSs;			// Vertex solid cross-sectional area
    d vVf;			// Vertex cell fluid volume
    d vVs;			// Vertex cell solid volume

    d SfR;			// Cross-sectional area of right face
    d SfL;			// Cross-sectional area of left  face

    d xR;			// Position of right cell wall
    d xL;			// Position of left cell wall
    d dxp;			// Distance to nearby right node
    d dxm;			// Distance to nearby left node
    d wLl,wRr,wLr,wRl;		// Weight functions for equations
    d wL0,wL1,wRNm1,wRNm2;    	// Special boundary weight functions

    
    virtual vd Error()=0;
    virtual dmat Jac();		// Returns the local Jacobian of this equation
    dmat zero;			// Zeros matrix of right size
    virtual dmat drhoim2();	// Derivative of current equation to density at node i-2
    virtual dmat dUim2();	// Etc
    virtual dmat dTim2();
    virtual dmat dpim2();
    virtual dmat dTsim2();

    
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

    virtual dmat drhoip2();	// Derivative of current equation to pressure
    virtual dmat dUip2();
    virtual dmat dTip2();
    virtual dmat dpip2();
    virtual dmat dTsip2();

    dmat D_r(); 		// Artificial viscosity right side
    dmat D_l();			// Artificial viscosity left size
    vd nu();			// Function of d^2p/dx^2
  };				// class Equation

} // namespace tube

