/*
 * duct.h
 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#pragma once
#ifndef DUCT_H_
#define DUCT_H_
#include "seg.h"
#include "constants.h"
#include "drag.h"
#include "heat.h"

namespace tasystem{
  class Jacobian;
  class TaSystem;
}

namespace duct{
  #ifndef SWIG
  SPOILNAMESPACE
  class DragResistance;
  class HeatSource;
  class Cell;
  class Geom;
  class BcCell;
  #endif  // SWIG
  #ifdef SWIG
  %catches(std::exception,...) Duct::setPhaseContraint(tasystem::PhaseConstraint v);
  %catches(std::exception,...) Duct::getValue(Varnr,us freqnr) const;
  %catches(std::exception,...) Duct::getValueC(Varnr,us freqnr) const;
  %catches(std::exception,...) Duct::getValueT(Varnr,d) const;
  %catches(std::exception,...) Duct::getErrorAt(us eqnr,us freqnr) const;
  %catches(std::exception,...) Duct::getCell(int i) const;
  #endif
  class Duct:public segment::Seg {
    // Pointer to possible PhaseConstraint instance
    tasystem::PhaseConstraint* pc_=nullptr;
    void showVertices(us detailnr) const ;   
    // Pointer to the geometry
    Geom* geom_=nullptr;		
  protected:
    std::vector<Cell*> cells;

  protected:
    // Copy constructor
    Duct(const Duct& other)=delete;
    Duct(const Duct& other,const tasystem::TaSystem& sys);
    Duct(const Geom& geom);
    Duct& operator=(const Duct&)=delete; // no copies allowed
  public:
    virtual ~Duct();          // Define this class as abstract

    void init();		// Initialize Duct

    // Push the right variables and equations
    virtual void setVarsEqs(Cell&) const;
    
    // Return a reference to the Geometry instance
    const Geom& geom() const;

    // Set individual at certain location for certain harmonic number
    void setResVar(Varnr,us i,us freqnr,d value);
    void setResVar(Varnr,us freqnr,const vd& value);
    // Optional: set a phase constraint on this segment,
    // somewere. Throws when wrong type of variable is constrained, or
    // when the freq number is not OK.
    void setPhaseContraint(tasystem::PhaseConstraint v);
    // Return a vector of all vertex positions
    vd getx() const;

    // Extract a result vector for given variable number
    // (rho,U,T,p,Ts) and frequency number.
    vd getValue(Varnr,us freqnr) const;

    // get result at certain time instance, normalized with period
    vd getValueT(Varnr,d timeinst) const;

    // Extract a  result vector  for given variable  number (rho,U,T,p,Ts)
    // and frequency number.
    vc getValueC(Varnr,us freqnr) const;

    // Extract a result vector for given variable number
    // (rho,U,T,p,Tw) and frequency number.
    vd getErrorAt(us eqnr,us freqnr) const;

    // Return number of Cells 
    us getNCells() const;

    // One way or another, Swig does not inherit this method to the
    // interface of a Duct. Therefore, we wrap it manually. Probably,
    // this should be improved later on.
    const char* __repr__() const {return segment::SegConBase::__repr__();}

    // Computes the Dof to constrain when called by a EngineSystem
    int providePhaseDof() const;
    d phaseDofValue() const;

    virtual vd error() const;

    // Methods not exposed to swig
    #ifndef SWIG
    void show(us showvertices=0) const;

    
    // Return number of DOFS
    us getNDofs() const;
    // Return # equations
    us getNEqs() const;    
    // Set dof nrs
    void setDofNrs(us firstdofnr);
    void setEqNrs(us firstdofnr);    
	// Obtain current mass in system
    d getMass() const;
    // Derivative of current mass in system to all dofs.
    void dmtotdx(vd&) const;
    // Set all higher harmonic amplitudes to zero
    void resetHarmonics();
    // Derivative of dofs to frequency
    void domg(vd& tofill) const;
    // Obtain result vector
    vd getRes() const;
    // Fill Jacobian submatrix
    void jac(tasystem::Jacobian& tofill) const;
    // Set result vector
    void setRes(const vd& res);
    // Update number of frequencies (something changed in Gc)
    void updateNf();    
    // ******************** End overloaded virtuals

    // return drag coefficient
    // virtual vd dragCoefVec(us) const=0;

    // Get reference to one of the boundary Cells
    const BcCell& bcCell(Pos) const;
    #endif
    // Access to cells

    // Can handle negative numbers
    // (-1 is last Cell
    const Cell& getCell(int i) const;
    #ifndef SWIG
    // If Duct has solid. Wall temperature is one of the variables
    // being solved. Otherwise not.
    virtual bool hasSolid() const {return false;}
    const Cell& operator[](us i) const;
    virtual const drag::DragResistance& dragResistance() const=0;
    virtual const HeatSource& heatSource() const=0;
    #endif
  private:
    void cleanup_cells();
  };				// Duct class

  #ifndef SWIG  
  inline Duct& asDuct(segment::Seg& s){return dynamic_cast<Duct&>(s);}
  inline const Duct& asDuct_const(const segment::Seg& s){return dynamic_cast<const Duct&>(s);}
  #endif
} /* namespace duct */

#endif /* DUCT_H_ */






