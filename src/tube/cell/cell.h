// File cell.h
#pragma once
#ifndef _CELL_H_
#define _CELL_H_

#include "var.h"
#include "tubeequation.h"
#include <map>
#include "constants.h"

namespace tasystem{class Jacobian;}
namespace tube {    

  #ifndef SWIG
  SPOILNAMESPACE;
  class Tube;
  class Equation;
  #endif

  #ifdef  SWIG
  %nodefaultctor Cell;
  #endif //  SWIG
  class Cell {
    //Gridpoint at a position
    //in a Tube
    const Tube* tube=nullptr;
  public:
    const tasystem::Globalconf* gc=nullptr;
  protected:
    const Cell* left_=nullptr;
    const Cell* right_=nullptr;

    tasystem::var rho_;         // Density at vertex
    tasystem::var ml_;		// Mass flow at left cell wall
    // tasystem::var mHl_;		// Total enthalpy flow at left cell wall
    tasystem::var mu_;      // Momentum flow at vertex
    tasystem::var T_;		// Temperature
    tasystem::var p_;      // Pressure at left cell wall
    tasystem::var Ts_;		// Solid temperature

    vector<tasystem::var*> vars;
    std::map<EqType,Equation*> eqs; // Vector of pointers to the equations

  public:
    // Geometric data ********************
    us i=0;                       // number of this cell
    d vx=0;                       // Vertex position
    d xl=0;                     // Absolute position of left cell wall
    d xr=0;                     // Absolute position of right cell wall    

    d vSf=0;			// Cell fluid cross-sectional area
    d vSs=0;			// Cell solid cross-sectional area
    d vVf=0;			// Cell cell fluid volume
    d vVs=0;			// Cell cell solid volume

    d Sfl=0,Sfr=0;		// Fluid surface area at cell walls.

    d Ssl=0,Ssr=0;    
    d vrh=0;			// Current cell hydraulic radius
    d rhl=0;            // Hydraulic radius at left cell wall
    d rhr=0;            // Hydraulic radius at right cell wall
    // End geometric data ********************

    #ifndef SWIG
    Cell(us i,const Tube&);
    virtual ~Cell();

    // No copy assignments allowed
    Cell& operator=(const Cell& v)=delete;

    // No copy constructors
    Cell(const Cell& )=delete;


    virtual void init(const Cell* left,const Cell* right);   

    const Cell* left() const {return left_;}
    const Cell* right() const {return right_;}

    const Tube& getTube() const {return *tube;}
    us geti() const {return i;}
    Equation* Eq(EqType et) {return eqs.at(et);}
    // Momentum flow at vertex position
    const tasystem::var& mu() const {return mu_;}
    const tasystem::var& ml() const {return ml_;}
    // const tasystem::var& mHl() const {return mHl_;}
    virtual const tasystem::var& mr() const {assert(right_); return right_->ml();}
    // virtual const tasystem::var& mHr() const { assert(right_); return right_->mHl();}
    const tasystem::var& p() const{return p_;}

    // These are variables for the left and right vertices, but are on
    // the cell walls for the leftmost and rightmost vertices
    // respectively
    const tasystem::var& rho() const {return rho_;}
    // virtual const tasystem::var& rhoL() const {assert(left_); return left_->rho();}
    // virtual const tasystem::var& rhoR() const {assert(right_); return right_->rho();}

    // const tasystem::var& U() const {return U_;}
    // virtual const tasystem::var& UL() const {assert(left_); return left_->U();}
    // virtual const tasystem::var& UR() const {assert(right_); return right_->U();}

    const tasystem::var& T() const { return T_;}
    virtual const tasystem::var& TL() const { assert(left_); return left_->T();}
    virtual const tasystem::var& TR() const {assert(right_); return right_->T();}
    const tasystem::var& pL() const { assert(left_); return left_->p();}
    const tasystem::var& pR() const {assert(right_); return right_->p();}
    const tasystem::var& rhoL() const { assert(left_); return left_->rho();}
    const tasystem::var& rhoR() const {assert(right_); return right_->rho();}

    const tasystem::var& Ts() const {return Ts_;}
    virtual const tasystem::var& TsL() const {assert(left_); return left_->Ts();}
    virtual const tasystem::var& TsR() const {assert(right_); return right_->Ts();}

    // Set this cell to be isentropically
    void setIsentropic();

    // Resets all higher harmonics. Can throw
    void resetHarmonics();

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
    void setResVar(Varnr,const tasystem::var& res);
    virtual void setResVar(Varnr,const vd& res); // Overridden for
                                                 // lefttubecell and
                                                 // righttubecell!
    
    virtual vd getRes() const;			  // Extract current result
                                          // vector
    d getValue(Varnr,us freqnr) const;

    #endif
    // Exposed to SWIG
    tasystem::var getValue(Varnr) const;
    #ifndef SWIG
    virtual void updateNf();

    // Convenience function, we need a lot of static (background
    // pressure) addings in the equations.
    // d Htot() const { return e.Htot();}
    d getMass() const;

    // These virtual functions are required such that boundary
    // condition sources can be added in a later stage by inheriting
    // from this Cell. By default these sources are not a
    // function of the dependent variables. That is why we do not have
    // to add Jacobian terms.
    virtual vd csource() const;	// Continuity source
    virtual vd msource() const;	// Momentum source
    virtual vd esource() const;	// Energy source
    #endif
  };				// Cell class
  
} // namespace tube

#endif /* _CELL_H_ */


