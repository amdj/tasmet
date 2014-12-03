// File tubevertex.h
#pragma once
#ifndef _TUBEVERTEX_H_
#define _TUBEVERTEX_H_

#include "tubeequation.h"
#include "localgeom.h"

#include "continuityeq.h"
#include "momentumeq.h"
#include "energyeq.h"
#include "stateeq.h"
#include "solidenergyeq.h"
#include "isentropiceq.h"
#include "var.h"

#include "varnr.h"

namespace tube{    


  SPOILNAMESPACE

  class Tube;

  // Abstract base class Vertex contains:
  // i: vertex nr
  // gc: pointer to Globalconf

  class TubeVertex{ //Gridpoint at a position
    //in a Tube
    us i;                       // number of this vertex
    LocalGeom lg;
    const Tube* tube=NULL;
    const tasystem::Globalconf* gc=NULL;
    const TubeVertex* left_=NULL;
    const TubeVertex* right_=NULL;

    vector<variable::var*> vars;
    vector<TubeEquation*> eqs; // Vector of pointers to the
  public:
    variable::var rho;		// Density
    variable::var U;		// Volume flow
    variable::var T;		// Temperature
    variable::var p;      // Pressure at left cell wall
    variable::var Ts;		// Solid temperature
  protected:
    // equations to solve for.
    Continuity c;
    Momentum m;
    Energy e;
    StateL sL;
    State s;
    SolidTPrescribed se;
    Isentropic is;              // Do we really need this burden?


    TubeVertex& operator=(const TubeVertex& v); // No copy assignments
    // allowed
    TubeVertex(const TubeVertex& );             // No copy constructors
  public:

    TubeVertex(us i,const Tube&);
    virtual ~TubeVertex(){}     // Will be derived from
    virtual void init(const TubeVertex* left,const TubeVertex* right);   

    // Get methods
    const LocalGeom& localGeom() const {return lg;}
    virtual const variable::var& pL() const;
    virtual const variable::var& pR() const;
    const TubeVertex* left() const {return left_;}
    const TubeVertex* right() const {return right_;}
    const Tube& getTube() const {return *tube;}
    us geti() const {return i;}

    void setIsentropic();
    void resetHarmonics();
    void setDofNrs(us firstdofnr);
    void setEqNrs(us firsteqnr);    
    us getNDofs() const;
    us getNEqs() const;

    virtual void setRes(const vd& res);			  // Set result vector
                                                  // to res
    // const methods
    virtual void show() const;
    virtual vd error() const;		       // Compute error for this gridpoint
    virtual void jac(Jacobian& tofill) const;		       // Fill complete Jacobian for this node
    virtual void domg(vd& ) const;
    void setRes(varnr,const variable::var& res);
    virtual vd getRes() const;			  // Extract current result
                                          // vector
    d getRes(varnr,us freqnr) const;
    variable::var getRes(varnr) const;
    virtual void updateNf();

    // Convenience function, we need a lot of static (background
    // pressure) addings in the equations.
    vd getp0t() const;

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


