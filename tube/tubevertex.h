// File tubevertex.h
#pragma once
#ifndef _TUBEVERTEX_H_
#define _TUBEVERTEX_H_

#include "vertex.h"
#include "tubeequation.h"
#include "localgeom.h"
#include "w.h"
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
  public:
    W::W w;
    dmat zero;			// Zeros matrix of right size
    const TubeVertex* left=NULL;
    const TubeVertex* right=NULL;
    us nCells;    

    // Equation-specific weight factors
    d cWddt=0,cWim1=0,cWi=0,cWip1=0;
    // Continuity artificial viscosity weight factor
    d cWart1=0,cWart2=0,cWart3=0,cWart4=0;		       
    d mWart1=0,mWart2=0,mWart3=0,mWart4=0;		       
    d mWddt=0,mWuim1=0,mWui=0,mWuip1=0;
    d mWpL=0,mWpR=0;
    d eWddt=0,eWddtkin=0,eWgim1=0,eWgi=0,eWgip1=0,eWc1=0,eWc2=0,eWc3=0,eWc4=0;      
    d eWkini=0,eWkinim1=0,eWkinip1=0;

    d eWispL=0,eWispR=0,eWisrho=0;
    
    // This is also the order in which they appear in the variable ptr
    // vector.
    variable::var rho;		// Density
    variable::var U;		// Volume flow
    variable::var T;		// Temperature
    variable::var p;      // Reference to p
    variable::var Ts;		// Solid temperature
    
    virtual const variable::var& pL() const;
    virtual const variable::var& pR() const;    

    std::vector<variable::var*> vars;

    vector<std::unique_ptr<TubeEquation> > eqs; // Vector of pointers to the
    // equations to solve for.

    virtual us getNDofs() const;
    virtual us getNEqs() const;
    void setDofNrs(us firstdofnr);
    void setEqNrs(us firstdofnr);    
    virtual ~TubeVertex(){}
    TubeVertex(){}
    TubeVertex(const TubeVertex& other){}
    virtual void setLeft(const Vertex&);
    virtual void setRight(const Vertex&);
    virtual void show() const;
    virtual vd error() const;		       // Compute error for this gridpoint
    virtual void jac(Jacobian& tofill) const;		       // Fill complete Jacobian for this node
    virtual void setRes(vd res);			  // Set result vector to res
    virtual vd domg() const;
    virtual vd getRes() const;			  // Extract current result vector
    // Convenience function, we need a lot of static (background
    // pressure) addings in the equations.
    vd getp0t() const;
  private:
    void connectTubeLeft(const SegBase& thisseg);
    void connectTubeRight(const SegBase& thisseg);
    void updateW(const SegBase& thisseg);
    void updateWEqs(const SegBase& thisseg);
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

    
  };				// TubeVertex class
} // namespace tube

#endif /* _TUBEVERTEX_H_ */


