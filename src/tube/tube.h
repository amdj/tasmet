/*
 * tube.h
 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#pragma once
#ifndef TUBE_H_
#define TUBE_H_
#include "seg.h"
#include "constants.h"
#include "drag.h"
#include "heat.h"

namespace tasystem{
  class Jacobian;
  class TaSystem;
}

namespace tube{
  #ifndef SWIG
  SPOILNAMESPACE
  class DragResistance;
  class HeatSource;
  class Cell;
  class Geom;
  class BcCell;
  #endif  // SWIG
  #ifdef SWIG
  %catches(std::exception,...) Tube::setPhaseContraint(tasystem::PhaseConstraint v);
  %catches(std::exception,...) Tube::getValue(Varnr,us freqnr) const;
  %catches(std::exception,...) Tube::getValueC(Varnr,us freqnr) const;
  %catches(std::exception,...) Tube::getErrorAt(us eqnr,us freqnr) const;
  #endif
  class Tube:public segment::Seg {
    // Pointer to possible PhaseConstraint instance
    tasystem::PhaseConstraint* pc_=nullptr;
    void showVertices(us detailnr) const ;   
    // Pointer to the geometry
    Geom* geom_=nullptr;		
  protected:
    std::vector<Cell*> cells;

  protected:
    // Copy constructor
    Tube(const Tube& other)=delete;
    Tube(const Tube& other,const tasystem::TaSystem& sys);
    Tube(const Geom& geom);
    Tube& operator=(const Tube&)=delete; // no copies allowed
  public:
    virtual ~Tube();          // Define this class as abstract

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

    // Extract a  result vector  for given variable  number (rho,U,T,p,Ts)
    // and frequency number.
    vc getValueC(Varnr,us freqnr) const;

    // Extract a result vector for given variable number
    // (rho,U,T,p,Ts) and frequency number.
    vd getErrorAt(us eqnr,us freqnr) const;

    // Return number of Cells 
    us getNCells() const;

    // One way or another, Swig does not inherit this method to the
    // interface of a Tube. Therefore, we wrap it manually. Probably,
    // this should be improved later on.
    const char* __repr__() const {return segment::SegConBase::__repr__();}

    // Computes the Dof to constrain when called by a EngineSystem
    int providePhaseDof() const;
    d phaseDofValue() const;

    // Methods not exposed to swig
    #ifndef SWIG

    virtual vd error() const;
    void show(us showvertices=0) const;


    us getNCell() const {return cells.size();}    

    // ******************** Overloaded virtual methods
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
    // Access to cells
    const Cell& getCell(us i) const;
    const Cell& operator[](us i) const;
    
    virtual const drag::DragResistance& getDragResistance() const=0;
    virtual const HeatSource& getHeatSource() const=0;
    #endif
  private:
    void cleanup_cells();
  };				// Tube class

  #ifndef SWIG  
  inline Tube& asTube(segment::Seg& s){return dynamic_cast<Tube&>(s);}
  inline const Tube& asTube_const(const segment::Seg& s){return dynamic_cast<const Tube&>(s);}
  #endif
} /* namespace tube */

#endif /* TUBE_H_ */






