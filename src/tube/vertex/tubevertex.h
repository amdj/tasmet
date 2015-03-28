// File tubevertex.h
#pragma once
#ifndef _TUBEVERTEX_H_
#define _TUBEVERTEX_H_

#include "localgeom.h"
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
  // Abstract base class Vertex contains:
  // i: vertex nr
  // gc: pointer to Globalconf

  class TubeVertex {
    //Gridpoint at a position
    //in a Tube
    us i;                       // number of this vertex
    const Tube* tube=NULL;
    WeightFactors* w_=NULL;
  public:
    const tasystem::Globalconf* gc=NULL;
  protected:
    const TubeVertex* left_=NULL;
    const TubeVertex* right_=NULL;

    vector<variable::var*> vars;
    vector<TubeEquation*> eqs; // Vector of pointers to the

    variable::var rho_;		// Density
    variable::var U_;		// Volume flow
    variable::var T_;		// Temperature
    variable::var pL_;      // Pressure at left cell wall
    variable::var Ts_;		// Solid temperature

  protected:
    // equations to solve for.
    Continuity c;
    Momentum m;
    Energy e;
    State s;
    SolidTPrescribed se;
    Isentropic is;              // Do we really need this burden?

  public:
    TubeVertex(us i,const Tube&);
    virtual ~TubeVertex();     // Deletes weightfactors instance

    // No copy assignments allowed
    TubeVertex& operator=(const TubeVertex& v)=delete;

    // No copy constructors
    TubeVertex(const TubeVertex& )=delete;

    const Continuity& continuity() const {return c;}
    const Momentum& momentum() const {return m;}
    const Energy& energy() const {return e;}


    virtual void init(const TubeVertex* left,const TubeVertex* right);   

    // Get methods
    const LocalGeom& localGeom() const;
    const WeightFactors& weightFactors() const;

    const TubeVertex* left() const {return left_;}
    const TubeVertex* right() const {return right_;}
    const Tube& getTube() const {return *tube;}
    us geti() const {return i;}

    const variable::var& pL() const { return pL_;}
    virtual const variable::var& pR() const {assert(right_); return right_->pL();}
    const variable::var p() const{return 0.5*(pL()+pR());}

    // These are variables for the left and right vertices, but are on
    // the cell walls for the leftmost and rightmost vertices
    // respectively
    const variable::var& rho() const {return rho_;}
    virtual const variable::var& rhoL() const {assert(left_); return left_->rho();}
    virtual const variable::var& rhoR() const {assert(right_); return right_->rho();}

    const variable::var& U() const {return U_;}
    virtual const variable::var& UL() const {assert(left_); return left_->U();}
    virtual const variable::var& UR() const {assert(right_); return right_->U();}

    const variable::var& T() const { return T_;}
    virtual const variable::var& TL() const { assert(left_); return left_->T();}
    virtual const variable::var& TR() const {assert(right_); return right_->T();}

    const variable::var& Ts() const {return Ts_;}
    virtual const variable::var& TsL() const {assert(left_); return left_->Ts();}
    virtual const variable::var& TsR() const {assert(right_); return right_->Ts();}

    // Set this vertex to be isentropically
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
    void setResVar(varnr,const variable::var& res);
    virtual void setResVar(varnr,const vd& res); // Overridden for
                                                 // lefttubevertex and
                                                 // righttubevertex!
    
    virtual vd getRes() const;			  // Extract current result
                                          // vector
    d getValue(varnr,us freqnr) const;
    variable::var getValue(varnr) const;
    virtual void updateNf();

    // Convenience function, we need a lot of static (background
    // pressure) addings in the equations.
    d Htot() const { return e.Htot();}
    d getCurrentMass() const;

    // These virtual functions are required such that boundary
    // condition sources can be added in a later stage by inheriting
    // from this TubeVertex. By default these sources are not a
    // function of the dependent variables. That is why we do not have
    // to add Jacobian terms.
    virtual vd csource() const;	// Continuity source
    virtual vd msource() const;	// Momentum source
    virtual vd esource() const;	// Energy source
  };				// TubeVertex class
} // namespace tube

#endif /* _TUBEVERTEX_H_ */


