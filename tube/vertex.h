// File vertex.h
#pragma once
#include "../var/var.h"
#include "geom.h"
#include "continuityeq.h"

#include "momentumeq.h"
#include "energyeq.h"
#include "stateeq.h"
#include "solidenergyeq.h"

#define Neq 5


namespace tube{
  
  class Vertex{
  public:
    Vertex(us i);
    virtual ~Vertex();
    const us i;			// The node number of this vertex
    
  };
  
    

  class TubeVertex:public Vertex{ //Gridpoint at a position in a Tube
  public:
    TubeVertex(const Tube& tube1,us i);
    ~TubeVertex();
    TubeVertex(const TubeVertex&); // The copy constructor.
    // TubeVertex* left;		// Left node pointer
    // TubeVertex* right;		// Right node pointer
 
    vd Error();				  // Compute error for this gridpoint
    dmat Jacobian();			  // Fill complete Jacobian for this node
    void Set(vd res);			  // Set result vector to res
    vd Get();				  // Extract current result vector

    const Tube& tube;			// Pointer to parent tube

    us Ns;
    d vSf;			// Vertex fluid cross-sectional area
    d vSs;			// Vertex solid cross-sectional area
    d vVf;			// Vertex cell fluid volume
    d vVs;			// Vertex cell solid volume
    d xR;			// Position of right cell wall
    d xL;			// Position of left cell wall
    d wLl,wRr,wLr,wRl;		// Weight functions for equations
    d wL0,wL1,wRNm1,wRNm2    	// Special boundary weight functions
    variable::var rho;		// Density
    variable::var U;		// Volume flow
    variable::var T;		// Temperature
    variable::var p;		// Pressure
    variable::var Ts;		// Solid temperature

    Equation* eq[Neq];		// Pointer array of all equations
    variable::var* vars[Neq];

    Continuity c;		// Continuity equation
    Momentum m;			// Momentum equation
    Energy e;			// Energy equation
    State s;			// State equation (ideal gas)
    Solidenergy se;		// Solid energy equation
  };				// TubeVertex class
} // namespace tube
