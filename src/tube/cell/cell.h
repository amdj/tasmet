// File cell.h
#pragma once
#ifndef _CELL_H_
#define _CELL_H_

#include "var.h"
#include "tubeequation.h"
#include "continuity.h"
#include "momentum.h"
#include "energy.h"
#include "state.h"
#include "solidenergy.h"
#include "isentropic.h"

#include "constants.h"

namespace tasystem{class Jacobian;}
namespace tube {    

  SPOILNAMESPACE;
  class Tube;
  // Abstract base class Cell contains:
  // i: cell nr
  // gc: pointer to Globalconf

  class Cell {
    //Gridpoint at a position
    //in a Tube
    const Tube* tube=nullptr;
  public:
    const tasystem::Globalconf* gc=nullptr;
  protected:
    const Cell* left_=nullptr;
    const Cell* right_=nullptr;

    vector<variable::var*> vars;
    vector<Equation*> eqs; // Vector of pointers to the

    variable::var rho_;		// Density
    variable::var rhoUL_;		// Mass flow at left cell wall
    variable::var T_;		// Temperature
    variable::var p_;      // Pressure at left cell wall
    variable::var Ts_;		// Solid temperature

  protected:
    // equations to solve for.
    Continuity c;
    Momentum m;
    // Temporarily put to isentropic
    // Energy e;
    Isentropic e;              // Do we really need this burden?
    State s;
    SolidTPrescribed se;
    Isentropic is;              // Does not have 

  public:
    // Geometric data ********************
    us i=0;                       // number of this cell
    d vx=0;                       // Vertex position
    d xL=0;                     // Absolute position of left cell wall
    d xR=0;                     // Absolute position of right cell wall    

    d vSf=0;			// Cell fluid cross-sectional area
    d vSs=0;			// Cell solid cross-sectional area
    d vVf=0;			// Cell cell fluid volume
    d vVs=0;			// Cell cell solid volume

    d SfL=0,SfR=0;		// Fluid surface area at cell walls.

    d SsL=0,SsR=0;    
    d vrh=0;			// Current cell hydraulic radius
    d rhL=0;            // Hydraulic radius at left cell wall
    d rhR=0;            // Hydraulic radius at right cell wall
    // End geometric data ********************

    Cell(us i,const Tube&);
    virtual ~Cell();

    // No copy assignments allowed
    Cell& operator=(const Cell& v)=delete;

    // No copy constructors
    Cell(const Cell& )=delete;

    const Continuity& continuity() const {return c;}
    const Momentum& momentum() const {return m;}
    // const Energy& energy() const {return e;}

    virtual void init(const Cell* left,const Cell* right);   

    const Cell* left() const {return left_;}
    const Cell* right() const {return right_;}
    const Tube& getTube() const {return *tube;}
    us geti() const {return i;}

    const variable::var& rhoUL() const {return rhoUL_;}
    vd U() const;
    virtual const variable::var& rhoUR() const {assert(right_); return right_->rhoUL();}
    const variable::var& p() const{return p_;}

    // These are variables for the left and right vertices, but are on
    // the cell walls for the leftmost and rightmost vertices
    // respectively
    const variable::var& rho() const {return rho_;}
    // virtual const variable::var& rhoL() const {assert(left_); return left_->rho();}
    // virtual const variable::var& rhoR() const {assert(right_); return right_->rho();}

    // const variable::var& U() const {return U_;}
    // virtual const variable::var& UL() const {assert(left_); return left_->U();}
    // virtual const variable::var& UR() const {assert(right_); return right_->U();}

    const variable::var& T() const { return T_;}
    virtual const variable::var& TL() const { assert(left_); return left_->T();}
    virtual const variable::var& TR() const {assert(right_); return right_->T();}

    const variable::var& Ts() const {return Ts_;}
    virtual const variable::var& TsL() const {assert(left_); return left_->Ts();}
    virtual const variable::var& TsR() const {assert(right_); return right_->Ts();}

    // Set this cell to be isentropically
    void setIsentropic();

    // Resets all higher harmonics. Can throw
    void resetHarmonics() throw (std::exception);

    void setDofNrs(us firstdofnr);
    void setEqNrs(us firsteqnr);    
    us getNDofs() const;
    us getNEqs() const;

    // Set result vector to res
    virtual void setRes(const vd& res);

    // Compute the error for all equations on this gridpoint
    vd error() const;
                                      

    // const methods
    virtual void show(us detailnr=1) const;
    vd errorAt(us i) const;

    virtual void jac(tasystem::Jacobian& tofill) const;		       // Fill complete Jacobian for this node
    virtual void domg(vd& ) const;
    void setResVar(Varnr,const variable::var& res);
    virtual void setResVar(Varnr,const vd& res); // Overridden for
                                                 // lefttubecell and
                                                 // righttubecell!
    
    virtual vd getRes() const;			  // Extract current result
                                          // vector
    d getValue(Varnr,us freqnr) const;
    variable::var getValue(Varnr) const;
    virtual void updateNf();

    // Convenience function, we need a lot of static (background
    // pressure) addings in the equations.
    // d Htot() const { return e.Htot();}
    d getCurrentMass() const;

    // These virtual functions are required such that boundary
    // condition sources can be added in a later stage by inheriting
    // from this Cell. By default these sources are not a
    // function of the dependent variables. That is why we do not have
    // to add Jacobian terms.
    virtual vd csource() const;	// Continuity source
    virtual vd msource() const;	// Momentum source
    virtual vd esource() const;	// Energy source
  };				// Cell class
} // namespace tube

#endif /* _CELL_H_ */


