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


namespace segment{
  
  class Vertex{
  public:
    Vertex(us i,const variable::varoperations& v);
    Vertex(const Vertex&);	// Copy constructor
    Vertex& operator=(const Vertex& v2); // Copy assignment
    virtual ~Vertex();
    
    const us i;			// The node number of this vertex
    const variable::varoperations& vop;
    const us& Ns;			// Number of sample points reference

    virtual vd Error();				  // Compute error for this gridpoint
    virtual dmat Jac();	       // Fill complete Jacobian for this node
    virtual void SetRes(vd res);			  // Set result vector to res
    virtual vd GetRes();				  // Extract current result vector

    tube::Equation* eq[Neq];		// Pointer array of all equations
    Vertex* left;
    Vertex* right;
    variable::var rho;		// Density
    variable::var U;		// Volume flow
    variable::var T;		// Temperature
    variable::var p;		// Pressure
    variable::var Ts;		// Solid temperature
    variable::var* vars[Neq];
  
  };
} // namespace segment

namespace tube{    

  class TubeVertex:public segment::Vertex{ //Gridpoint at a position in a Tube
  public:
    TubeVertex(const Tube& tube1,us i);
    TubeVertex(const TubeVertex&); // The copy constructor.
    void operator=(const TubeVertex&);
    virtual ~TubeVertex();

    
    // TubeVertex* left;		// Left node pointer
    // TubeVertex* right;		// Right node pointer
     
    const Tube& tube;			// Pointer to parent tube

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
    
  };				// TubeVertex class
} // namespace tube



