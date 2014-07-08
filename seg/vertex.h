// File vertex.h
#pragma once
#ifndef _VERTEX_H_
#define _VERTEX_H_

#include "var.h"
#include "globalconf.h"
#include "geom.h"
#include "equation.h"

#define Neq (5)


namespace segment{
  SPOILNAMESPACE
  using  tasystem::Globalconf;
  using variable::var;
  
  class Seg;  
  class Vertex{
  public:
    Vertex(const Seg& seg,us i);
    Vertex(const Vertex&);	// Copy constructor
    Vertex& operator=(const Vertex& v2); // Copy assignment
    void Init(const Globalconf& gc){
      this->gc=&gc;
      this->updateW();
      rho=var(gc);
      U=var(gc);
      T=var(gc);
      p=var(gc);
      Ts=var(gc);
    }
    virtual ~Vertex();
    const Seg& seg;
    const us i;			// The node number of this vertex
    const tasystem::Globalconf *gc=NULL;

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

    virtual void updateW();	       // Update weight functions of equations
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

    d vxim1,vxi,vxip1;    		// Vertex position of left and right vertex
  };
} // namespace segment


#endif /* _VERTEX_H_ */
