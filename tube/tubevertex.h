// File tubevertex.h
#pragma once
#ifndef _TUBEVERTEX_H_
#define _TUBEVERTEX_H_

#include "vertex.h"
#include "tubeequation.h"
#include "localgeom.h"

#include "continuityeq.h"
#include "momentumeq.h"
#include "energyeq.h"
#include "stateeq.h"
#include "solidenergyeq.h"
#include "isentropiceq.h"
#include "var.h"

namespace segment{class SegBase;}



namespace tube{    
  SPOILNAMESPACE
  using segment::SegBase;
  using segment::LocalGeom;
  class Tube;

  // Abstract base class Vertex contains:
  // i: vertex nr
  // gc: pointer to Globalconf
  // VertexVec left,right : vector of pointers to left and right
  // vertices. Reserved for later more complicated stuff.
  
  class TubeVertex:public segment::Vertex{ //Gridpoint at a position in a Tube
  protected: 
    d wLl=0,wRr=0,wLr=0,wRl=0;		// Basic weight functions
    d wL0=0,wL1=0,wRNm1=0,wRNm2=0;    	// Special boundary weight functions
    d xvim1=0,xvi=0,xvip1=0;
    d dxm=0,dxp=0;
    d vSfR=0,vSfL=0;		// Cross sectional area at x-position
    // of vertex.
    int UsignL=1;
    int UsignR=1;

  public:
    friend class StateR;
    friend class Energy;
    friend class RightImpedance;
    friend class RightImpedanceMomentumEq;
    friend class RightTwImpedanceEq;
    friend class LeftImpedance;
    friend class TwImpedance;    

    const TubeVertex* left=NULL;
    const TubeVertex* right=NULL;
    us nCells;    

    // This is also the order in which they appear in the variable ptr
    // vector.
    variable::var rho;		// Density
    variable::var U;		// Volume flow
    variable::var T;		// Temperature
    variable::var p;      // Reference to p
    variable::var Ts;		// Solid temperature
    
    virtual const variable::var& pL() const;
    virtual const variable::var& pR() const;

    vector<variable::var*> vars;
    vector<TubeEquation*> eqs; // Vector of pointers to the
    // equations to solve for.
    Continuity c;
    Momentum m;
    Energy e;
    StateL sL;
    State s;
    SolidTPrescribed se;
    Isentropic is;
    virtual void setpR(const variable::var& o) {
      WARN("pR tried to be set on normal tubeVertex!");
    }

    void setIsentropic();
    void resetHarmonics();
    virtual us getNDofs() const;
    virtual us getNEqs() const;
    void setDofNrs(us firstdofnr);
    void setEqNrs(us firstdofnr);    
    virtual ~TubeVertex(){}
    TubeVertex(){}
    TubeVertex(const TubeVertex& v){} // Copy does nothing
    TubeVertex& operator=(const TubeVertex& v){return *this;} // Assignment does nothing
    virtual void setLeft(const Vertex&);
    virtual void setRight(const Vertex&);
    virtual void show() const;
    virtual vd error() const;		       // Compute error for this gridpoint
    virtual void jac(Jacobian& tofill) const;		       // Fill complete Jacobian for this node
    virtual void setRes(vd res);			  // Set result vector to res
    virtual void domg(vd& ) const;
    virtual vd getRes() const;			  // Extract current result vector
    virtual void updateNf();
    // Convenience function, we need a lot of static (background
    // pressure) addings in the equations.
    vd getp0t() const;
  public:
    virtual void initTubeVertex(us i,const Tube&);   
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
   
    void connectTubeLeft(const Tube& thisseg);
    void connectTubeRight(const Tube& thisseg);
    void updateW(const Tube& thisseg);
    
  };				// TubeVertex class
} // namespace tube

#endif /* _TUBEVERTEX_H_ */


