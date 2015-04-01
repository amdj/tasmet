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
  #endif

  class Tube:public segment::Seg {

    void showVertices(us detailnr) const ;   
    // Pointer to the geometry
    Geom* geom_=nullptr;		
  protected:
    std::vector<Cell*> cells;

  protected:
    // Copy constructor
    Tube(const Tube& other);
    Tube(const Geom& geom);
    Tube& operator=(const Tube&)=delete; // no copies allowed
  public:
    virtual ~Tube();          // Define this class as abstract

    // Return a reference to the Geometry instance
    const Geom& geom() const;

    // Return the time-average total enthalpy flow (Watts)
    // vd Htot() const throw(std::exception);

    // Set individual at certain location for certain harmonic number
    void setResVar(Varnr,us i,us freqnr,d value);
    void setResVar(Varnr,us freqnr,const vd& value);
    
    // Return a vector of all dof positions, including the DOFS at the
    // ends.
    vd getx() const;
    
    // Return a value in form of an array with length equal to the
    // length returned with getx(). 
    vd getValue(Varnr,us freqnr) const throw(std::exception); // Extract a result vector for given variable number (rho,U,T,p,Ts) and frequency number.
    vc getValueC(Varnr,us freqnr) const throw(std::exception); // Extract a result vector for given variable number (rho,U,T,p,Ts) and frequency number.
    vd getErrorAt(us eqnr,us freqnr) const throw(std::exception); // Extract a result vector for given variable number (rho,U,T,p,Ts) and frequency number.

     us getNCells() const;
    // Methods not exposed to swig

    #ifndef SWIG
    virtual vd error() const;
    virtual void init(const tasystem::TaSystem&);
    void setRes(const segment::Seg& other); // To copy from a
    void show(us showvertices=0) const;


    us getNCell() const {return cells.size();}    
    vd interpolateResMid(Varnr v,d x) const; // Amplitude data vectors
    vd interpolateResStaggered(Varnr v,d x) const; // Amplitude data!!

    virtual us getNDofs() const;
    virtual us getNEqs() const;    
    virtual void setDofNrs(us firstdofnr);
    virtual void setEqNrs(us firstdofnr);    
    virtual d getCurrentMass() const;	// Obtain current mass in
                                        // system
    virtual void dmtotdx(vd&) const; // Derivative of current mass in
				    // system to all dofs.
    virtual void resetHarmonics();             // Set all higher
                                               // harmonic amplitudes
                                               // to zero
    virtual void domg(vd& tofill) const;

    virtual vd getRes() const;
    virtual d getRes(us dofnr) const;
    virtual void setRes(const vd& res);
    virtual void updateNf();    

    // *similar* segment
    virtual vd dragCoefVec(us) const=0;              // return drag
                                                   // coefficient

    const Cell& operator[](us i) const;
    const BcCell& bcCell(Pos) const;
    const Cell& getCell(us i) const;
    virtual void jac(tasystem::Jacobian& tofill) const;
    virtual const DragResistance& getDragResistance() const=0;
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






