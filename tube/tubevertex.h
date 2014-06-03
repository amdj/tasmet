// File vertex.h
#pragma once
#ifndef _TUBEVERTEX_H_
#define _TUBEVERTEX_H_

#include "vertex.h"
#include "continuityeq.h"
#include "momentumeq.h"
#include "energyeq.h"
#include "stateeq.h"
#include "solidenergyeq.h"



namespace tube{    
  class RightImpedanceMomentumEq;
  class LeftPressure;
  
  class TubeVertex:public segment::Vertex{ //Gridpoint at a position in a Tube
  public:
    TubeVertex(const Tube& tube1,us i);
    TubeVertex(const TubeVertex&); // The copy constructor.
    void operator=(const TubeVertex&);
    virtual ~TubeVertex();

    const Tube& tube;			// Pointer to parent tube
    virtual void updateW();

    Continuity c;		// Continuity equation
    Momentum m;			// Momentum equation
    Energy e;			// Energy equation
    State s;			// State equation (ideal gas)
    Solidenergy se;		// Solid energy equation
    Isentropic is;
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
    
  protected:
    d vSf;			// Vertex fluid cross-sectional area
    d vSs;			// Vertex solid cross-sectional area
    d vVf;			// Vertex cell fluid volume
    d vVs;			// Vertex cell solid volume

    d SfR;			// Cross-sectional area of right face
    d SfL;			// Cross-sectional area of left  face

    d xR;			// Position of right cell wall
    d xL;			// Position of left cell wall
    d dxp;			// Distance to nearby right node
    d dxm;			// Distance to nearby left node
    d wLl,wRr,wLr,wRl;		// Weight functions for equations
    d wL0,wL1,wRNm1,wRNm2;    	// Special boundary weight functions
    
  };				// TubeVertex class
} // namespace tube

#endif /* _TUBEVERTEX_H_ */


