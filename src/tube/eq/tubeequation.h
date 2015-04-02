// tubeequation.h
//
// Author: J.A. de Jong 
//
// Description:
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef TUBEEQUATION_H
#define TUBEEQUATION_H
#include "vtypes.h"
#include "constants.h"

namespace tasystem{
  class JacRow;
  class JacCol;
}
namespace variable{class var;}

namespace tube{
  
  SPOILNAMESPACE  
  class Cell;
  class Tube;
  class WeightFactors;
  
  class Equation{
  protected:
    const Cell& v;
    us dofnr;
  public:
    Equation(const Cell& v):v(v){TRACE(15,"Equation(v)");}
    void setDofNr(us Dofnr){dofnr=Dofnr;}
    us getDofNr(){return dofnr;}    
    virtual void init()=0;
    virtual enum EqType getType() const=0;

    // Return an eye of the right size:
    dmat eye() const;
    vd zeros() const;
    virtual tasystem::JacRow jac() const=0;		// Returns the local Jacobian of this equation
    virtual vd error() const=0;
    virtual void show() const { cout << "Empty equation description. From equation.h.\n";}
    virtual void domg(vd&) const {/* Placeholder */}
    // The definition of these factors is the original definition of d_j+/-1/2 of Wesseling:
    // Wesseling: d_j+/-1/2 = r_j+/-0.5 * epsilon_j+/-0.5
    // dmat d_r() const; 		// Artificial viscosity pre-factor right side 
    // dmat d_l() const;			// Artificial viscosity pre-factor
                                // left size
    vd getp0t() const;
    vd getT0t() const;
    virtual ~Equation(){}
  private:
    vd nu(const Cell&) const;			// Function of d^2p/dx^2
  };				// class Equation
} // namespace tube
#endif /* _TUBEEQUATION_H_ */

//////////////////////////////////////////////////////////////////////
