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

namespace segment{class SegBase;}



namespace tube{    


  SPOILNAMESPACE

  class Tube;

  // Abstract base class Vertex contains:
  // i: vertex nr
  // gc: pointer to Globalconf
  // VertexVec left,right : vector of pointers to left and right
  // vertices. Reserved for later more complicated stuff.
class WeightFunctions
  {
  public:
    d wLl=0,wRr=0,wLr=0,wRl=0;		// Basic weight functions
    d wL0=0,wL1=0,wRNm1=0,wRNm2=0;    	// Special boundary weight functions
    d vxm1=0,vx=0,vxp1=0;
    d dxm=0,dxp=0;
    d vSfR=0,vSfL=0;		// Cross sectional area at x-position
};  
  class TubeVertex{ //Gridpoint at a position
    //in a Tube
    const tasystem::Globalconf* gc=NULL;
    const Tube* tube;
    us i;
    LocalGeom lg;
    const TubeVertex* left=NULL;
    const TubeVertex* right=NULL;

    vector<variable::var*> vars;
    vector<TubeEquation*> eqs; // Vector of pointers to the

    WeightFunctions w_;
    // Later to be removed
    d wLl=0,wRr=0,wLr=0,wRl=0;		// Basic weight functions
    d wL0=0,wL1=0,wRNm1=0,wRNm2=0;    	// Special boundary weight functions
    d vxm1=0,vx=0,vxp1=0;
    d dxm=0,dxp=0;
    d vSfR=0,vSfL=0;		// Cross sectional area at x-position

    Continuity c;
    Momentum m;
    Energy e;
    StateL sL;
    State s;
    SolidTPrescribed se;
    Isentropic is;

    variable::var rho;		// Density
    variable::var U;		// Volume flow
    variable::var T;		// Temperature
    variable::var p;      // Pressure at left cell wall
    variable::var Ts;		// Solid temperature

    TubeVertex& operator=(const TubeVertex& v); // No copy assignments
                                                // allowed
    TubeVertex(const TubeVertex& );             // No copy constructors
  public:
    // This is also the order in which they appear in the variable ptr
    // vector.
    virtual const variable::var& pL() const;
    virtual const variable::var& pR() const;
    void setLeft(const TubeVertex&);
    void setRight(const TubeVertex&);

    // equations to solve for.

    void setIsentropic();
    void resetHarmonics();
    void setDofNrs(us firstdofnr);
    void setEqNrs(us firsteqnr);    
    us getNDofs() const;
    us getNEqs() const;
    TubeVertex(){}
    virtual ~TubeVertex(){}
    virtual void init(us i,const Tube&);   
    virtual void setRes(const vd& res);			  // Set result vector
                                                  // to res

    // const methods
    const WeightFunctions& w() const {return w_;}
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

    // These virtual functions are reqsuired such that boundary
    // condition sources can be added in a later stage by inheriting
    // from this TubeVertex. By default these sources are not a
    // function of the dependent variables. That is why we do not have
    // to add Jacobian terms.
    virtual vd csource() const;	// Continuity source
    virtual vd msource() const;	// Momentum source
    virtual vd esource() const;	// Energy source
  private:
    void leftVertex();
    void middleVertex();
    void rightVertex();
    void allVertex();
    void updateW(const Tube& thisseg);
    
  };				// TubeVertex class
} // namespace tube

#endif /* _TUBEVERTEX_H_ */


