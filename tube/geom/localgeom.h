#pragma once			// 
#ifndef _LOCALGEOM_H_
#define _LOCALGEOM_H_
#include "vtypes.h"

namespace tube {
  SPOILNAMESPACE

  class Geom;
  class LocalGeom
  {
  public:
    LocalGeom(const Geom& geom,us i);
    virtual ~LocalGeom(){}
    const Geom* geom=NULL;	// Pointer to global Geometry
    us i;			// Vertex number
    us nCells() const{return geom->nCells();}	

    d vSf=0;			// Vertex fluid cross-sectional area
    d vSs=0;			// Vertex solid cross-sectional area
    d vVf=0;			// Vertex cell fluid volume
    d vVs=0;			// Vertex cell solid volume

    d SfL=0,SfR=0;		// Fluid surface area at cell walls. *WARNING*
    // use these only for end boundary *conditions!!
    
    d vrh=0;			// Current vertex hydraulic radius
    d xR=0;			// Absolute position of right cell wall
    d xL=0;			// Absolute position of left cell wall
    d xr=0,xl=0;			// Postition of right and left wall
    // relative to vertex position
    
    d dxp=0;			// Distance to nearby right node
    d dxm;			// Distance to nearby left node

    d vx=0;
    // right vertex.

    virtual void show() const;

    // Geom can call constructor of this class
    friend class Geom;
  };				// class LocalGeom

  
} // namespace tube

#endif /* _LOCALGEOM_H_ */
