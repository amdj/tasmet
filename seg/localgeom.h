#pragma once			// 
#ifndef _LOCALGEOM_H_
#define _LOCALGEOM_H_
#include <vtypes.h>

namespace segment {
  SPOILNAMESPACE

  class Geom;
  class LocalGeom
  {
  private:
    LocalGeom(const Geom& geom,us i);
  public:
    LocalGeom(){ TRACE(10,"Created empty LocalGeom"); }
    virtual ~LocalGeom(){}
    us i;			// Vertex number
    us nCells;
    d vSf;			// Vertex fluid cross-sectional area
    d vSs;			// Vertex solid cross-sectional area
    d vVf;			// Vertex cell fluid volume
    d vVs;			// Vertex cell solid volume

    d SfR;			// Cross-sectional area of right face
    d SfL;			// Cross-sectional area of left  face

    d xR;			// Absolute position of right cell wall
    d xL;			// Absolute position of left cell wall
    d xr,xl;			// Postition of right and left wall relative to vertex position
    d dxp;			// Distance to nearby right node
    d dxm;			// Distance to nearby left node

    d vxim1,vxi,vxip1;    		// Vertex position of left and right vertex
    void show();
    friend class Geom;
  };				// class LocalGeom

  
}                               // namespace segment

#endif /* _LOCALGEOM_H_ */
