#pragma once			// 
#ifndef _LOCALGEOM_H_
#define _LOCALGEOM_H_
#include "vtypes.h"

namespace tube {

  #ifndef SWIG
  SPOILNAMESPACE
  #endif

  class Geom;
  class Cell;

  class LocalGeom
  {
  public:
    LocalGeom(const Cell&);
    LocalGeom(const Geom& geom,us i);
    virtual ~LocalGeom(){}
    const Geom* geom=NULL;	// Pointer to global Geometry
    us i;			// Cell number
    us nCells() const;

    d vSf=0;			// Cell fluid cross-sectional area
    d vSs=0;			// Cell solid cross-sectional area
    d vVf=0;			// Cell cell fluid volume
    d vVs=0;			// Cell cell solid volume

    d SfL=0,SfR=0;		// Fluid surface area at cell walls. *WARNING*
    // use these only for end boundary *conditions!!
    d SsL=0,SsR=0;    
    d vrh=0;			// Current cell hydraulic radius
    d xR=0;			// Absolute position of right cell wall
    d xL=0;			// Absolute position of left cell wall

    d vx=0;
    // right cell.

    virtual void show() const;

    // Geom can call constructor of this class
    friend class Geom;
  };				// class LocalGeom

  
} // namespace tube

#endif /* _LOCALGEOM_H_ */
