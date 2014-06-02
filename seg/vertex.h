// File vertex.h
#pragma once
#ifndef _VERTEX_H_
#define _VERTEX_H_

#include "var.h"
#include "geom.h"
#include "equation.h"



#define Neq (5)


namespace segment{
  SPOILNAMESPACE
  
  class Seg;  
  class Vertex{
  public:
    Vertex(const Seg& seg,us i);
    Vertex(const Vertex&);	// Copy constructor
    Vertex& operator=(const Vertex& v2); // Copy assignment
    virtual ~Vertex();
    const Seg& seg;
    const us i;			// The node number of this vertex
    const tasystem::Globalconf& gc;
    const us& Ns;			// Number of sample points reference

    virtual vd Error();				  // Compute error for this gridpoint
    virtual dmat Jac();	       // Fill complete Jacobian for this node
    virtual void SetRes(vd res);			  // Set result vector to res
    virtual vd GetRes();				  // Extract current result vector

    segment::Equation* eq[Neq];		// Pointer array of all equations
    const Vertex* left;
    const Vertex* right;
    variable::var rho;		// Density
    variable::var U;		// Volume flow
    variable::var T;		// Temperature
    variable::var p;		// Pressure
    variable::var Ts;		// Solid temperature
    variable::var* vars[Neq];
  
  };
} // namespace segment


#endif /* _VERTEX_H_ */
