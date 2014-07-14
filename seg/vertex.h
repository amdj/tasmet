// File vertex.h
#pragma once
#ifndef _VERTEX_H_
#define _VERTEX_H_

#include "segbase.h"
#include "var.h"
#include "equation.h"

#define Neq (5)


namespace segment{
  SPOILNAMESPACE
  using tasystem::Globalconf;
  using segment::Geom;
  using variable::var;

  
  class Vertex{
  public:
    us i=0;			// The node number of this vertex
    us Ncells=0;
    const Globalconf *gc=NULL;
    Equation* eq[Neq];		// Pointer array of all equations

    const Vertex* left=NULL;
    const Vertex* right=NULL;
    variable::var rho;		// Density
    variable::var U;		// Volume flow
    variable::var T;		// Temperature
    variable::var p;		// Pressure
    variable::var Ts;		// Solid temperature
    variable::var* vars[Neq]={&rho,&U,&T,&p,&Ts};
    
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
    
    Vertex();
    virtual void show();
    virtual ~Vertex();
    // Standard copy constructor will suffice
    Vertex(const Vertex&); 
    Vertex& operator=(const Vertex& v2); // Copy assignment
    virtual void Init(us i,const Globalconf& gc,const Geom& geom);
  private:
    void updateW(const Geom&);	       // Update weight functions of equations
  public:
    vd Error();				  // Compute error for this gridpoint
    dmat Jac();	       // Fill complete Jacobian for this node
    void SetRes(vd res);			  // Set result vector to res
    vd GetRes();				  // Extract current result vector




  };
} // namespace segment


#endif /* _VERTEX_H_ */
