#pragma once
#ifndef _TUBEEQUATION_H_
#define _TUBEEQUATION_H_
#define Neq (5)
#include "equation.h"


namespace tube{
  SPOILNAMESPACE  
  using namespace segment;
  using tasystem::Globalconf;
  using segment::Equation;

  class TubeVertex;
  class TubeEquation:public Equation{
  public:
      // const dmat& fDFT,iDFT,DDTfd;	// forward, backward dicrete fourier transform, derivative to time matrix (freq domain)
    // ****************************** THIS METHOD SHOULD NOT BE OVERRIDDEN!!!
    dmat jac(const TubeVertex& v) const;		// Returns the local Jacobian of this equation
    // ******************************
    virtual vd error(const TubeVertex&) const=0;
    // The definition of these factors is the original definition of d_j+/-1/2 of Wesseling:
    // Wesseling: d_j+/-1/2 = r_j+/-0.5 * epsilon_j+/-0.5
    dmat d_r(const TubeVertex&) const; 		// Artificial viscosity pre-factor right side 
    dmat d_l(const TubeVertex&) const;			// Artificial viscosity pre-factor left size
  private:
    vd nu(const TubeVertex&) const;			// Function of d^2p/dx^2

  public:
    virtual dmat drhoim2(const TubeVertex&) const;	// Derivative of current equation to density at node i-2
    virtual dmat dUim2(const TubeVertex&) const;	// Etc
    virtual dmat dTim2(const TubeVertex&) const;
    virtual dmat dpim2(const TubeVertex&) const;
    virtual dmat dTsim2(const TubeVertex&) const;
    
    virtual dmat drhoim1(const TubeVertex&) const;	// Derivative of current equation to density at node i-1
    virtual dmat dUim1(const TubeVertex&) const;	// Etc
    virtual dmat dTim1(const TubeVertex&) const;
    virtual dmat dpim1(const TubeVertex&) const;
    virtual dmat dTsim1(const TubeVertex&) const;

    virtual dmat drhoi(const TubeVertex&) const;	// Derivative of current equation to density at node i
    virtual dmat dUi(const TubeVertex&) const;		// Etc
    virtual dmat dTi(const TubeVertex&) const;
    virtual dmat dpi(const TubeVertex&) const;
    virtual dmat dTsi(const TubeVertex&) const;

    virtual dmat drhoip1(const TubeVertex&) const;	// Derivative of current equation to pressure
    virtual dmat dUip1(const TubeVertex&) const;
    virtual dmat dTip1(const TubeVertex&) const;
    virtual dmat dpip1(const TubeVertex&) const;
    virtual dmat dTsip1(const TubeVertex&) const;

    virtual dmat drhoip2(const TubeVertex&) const;	// Derivative of current equation to pressure
    virtual dmat dUip2(const TubeVertex&) const;
    virtual dmat dTip2(const TubeVertex&) const;
    virtual dmat dpip2(const TubeVertex&) const;
    virtual dmat dTsip2(const TubeVertex&) const;
  };				// class TubeEquation
} // namespace tube
#endif /* _TUBEEQUATION_H_ */
