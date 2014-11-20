#pragma once			// 
#ifndef _LOCALGEOM_H_
#define _LOCALGEOM_H_
#include "vtypes.h"

namespace geom {
  SPOILNAMESPACE

  class Geom;
  class LocalGeom
  {
  private:
    LocalGeom(const Geom& geom,us i);
  public:
    LocalGeom(){ TRACE(10,"Created empty LocalGeom"); }

    const Geom* geom=NULL;	// Pointer to global Geometry
    us i;			// Vertex number
    us nCells;			// Total number of Cells in this segment
    
    d vSf;			// Vertex fluid cross-sectional area
    d vSs;			// Vertex solid cross-sectional area
    d vVf;			// Vertex cell fluid volume
    d vVs;			// Vertex cell solid volume

    d SfL,SfR;		// Fluid surface area at cell walls. *WARNING*
    // use these only for end boundary
    // *conditions!!
    
    d vrh;			// Current vertex hydraulic radius
    d xR;			// Absolute position of right cell wall
    d xL;			// Absolute position of left cell wall
    d xr,xl;			// Postition of right and left wall
				// relative to vertex position
    
    d dxp;			// Distance to nearby right node
    d dxm;			// Distance to nearby left node

    d vxi;
    // right vertex.

    void show() const;

    // Geom can call constructor of this class
    friend class Geom;
  };				// class LocalGeom

  
} // namespace geom

#endif /* _LOCALGEOM_H_ */
