// File tubeequation.h
#pragma once
#ifndef _TUBEEQUATION_H_
#define _TUBEEQUATION_H_
#include "vtypes.h"
#include "jacobian.h"

namespace tube{
  using namespace tasystem;
  
  enum EqType{
    Con,			// Continuity
    Mom,			// Momentum
    Ene,			// Energy
    Ise,			// Isentropic
    Sta,			// State
    Sol,			// SolidEnergy
    Non				// None
  };
  
  SPOILNAMESPACE  
  class Tube;
  class TubeVertex;
  class TubeEquation{
  protected:
    us dofnr;
  public:
    void setDofNr(us Dofnr){dofnr=Dofnr;}
    us getDofNr(){return dofnr;}    
    virtual void init(const Tube&) {}
    virtual enum EqType getType() const { return EqType::Non;}

    virtual JacRow jac(const TubeVertex& v) const=0;		// Returns the local Jacobian of this equation
    virtual vd error(const TubeVertex&) const=0;
    virtual void show() const { cout << "Empty equation description. From equation.h.\n";}
    virtual void domg(const TubeVertex&,vd&) const {/* Placeholder */}
    // The definition of these factors is the original definition of d_j+/-1/2 of Wesseling:
    // Wesseling: d_j+/-1/2 = r_j+/-0.5 * epsilon_j+/-0.5
    dmat d_r(const TubeVertex&) const; 		// Artificial viscosity pre-factor right side 
    dmat d_l(const TubeVertex&) const;			// Artificial viscosity pre-factor left size
    virtual ~TubeEquation(){}
  private:
    vd nu(const TubeVertex&) const;			// Function of d^2p/dx^2
  };				// class TubeEquation
} // namespace tube
#endif /* _TUBEEQUATION_H_ */
