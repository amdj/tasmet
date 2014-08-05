// File vertex.h
#pragma once
#ifndef _TUBEVERTEX_H_
#define _TUBEVERTEX_H_

#include "tube.h"
#include "vertex.h"



namespace tube{    
  SPOILNAMESPACE
  class RightImpedanceMomentumEq;
  class LeftPressure;
  using segment::Geom;
  class TubeEquation;
  using tasystem::Globalconf;
  
  class TubeVertex:public segment::Vertex{ //Gridpoint at a position in a Tube

  public:
    // protected:
    dmat zero;			// Zeros matrix of right size
    us i,nCells;
    const TubeVertex* left=NULL;
    const TubeVertex* right=NULL;
    d wLl,wRr,wLr,wRl;		// Weight functions for equations
    d wL0,wL1,wRNm1,wRNm2;    	// Special boundary weight functions
    const d K=10.0;
    d cWddt,cWim1,cWi,cWip1;
    d mWddt,mWuim1,mWui,mWuip1,mWpim1,mWpi,mWpip1;
    d eWddt,eWgim1,eWgi,eWgip1,eWjim1,eWji,eWjip1,eWc1,eWc2,eWc3,eWc4;      
    d eWkini,eWkinim1,eWkinip1;
    
    variable::var rho;		// Density
    variable::var U;		// Volume flow
    variable::var T;		// Temperature
    variable::var p;		// Pressure
    variable::var Ts;		// Solid temperature
    vector<variable::var*> vars={&rho,&U,&T,&p,&Ts};
    vector<const TubeEquation*> eq;

    TubeVertex();
    TubeVertex(const TubeVertex&);
    TubeVertex& operator=(const TubeVertex&);
    virtual ~TubeVertex();
    virtual void show() const;
    virtual vd error() const;		       // Compute error for this gridpoint
    virtual dmat jac() const;		       // Fill complete Jacobian for this node
    virtual void setRes(vd res);			  // Set result vector to res

    virtual vd getRes() const;			  // Extract current result vector
    vd getp0t() const;
  private:
    void updateW(const SegBase& thisseg);
  public:
    virtual void initTubeVertex(us i,const Tube&);   
    // These virtual functions are required such that boundary
    // condition sources can be added in a later stage by inheriting
    // from this TubeVertex. By default these sources are not a
    // function of the dependent variables. That is why we do not have
    // to add Jacobian terms.
    virtual vd csource() const;	// Continuity source
    virtual vd msource() const;	// Momentum source
    virtual vd esource() const;	// Energy source

    friend class TubeEquation;   
    friend class Continuity;
    friend class Momentum;
    friend class Energy;
    friend class State;    
    friend class Solidenergy;    
    friend class RightImpedanceMomentumEq;
    friend class LeftPressure;
    
    
  };				// TubeVertex class
} // namespace tube

#endif /* _TUBEVERTEX_H_ */


