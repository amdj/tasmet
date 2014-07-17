#pragma once			// 
#ifndef _GEOM_H_
#define _GEOM_H_
#include <vtypes.h>
#define MAXGP 500

namespace segment {
  SPOILNAMESPACE
  class Geom{
  public:
    Geom(us gp,d L,d S,d phi,d rh,string cshape);
    ~Geom();
    const string shape;	 // Shape keyword: currently available: 'circ','vert','blapprox'
    void show();
    const us gp;         // Number of cell walls
    const us Ncells;	 // Number of cells
    d L;		 //Length of the Segment
    vd x;		 // Position of cell walls

    vd S;		 // Cross sectional area as a function of x
    vd Ss;		 // Solid-occupied cross-sectional area
    vd Sf;		 // Fluid-occupied cross-sectional area of cell-wall

    vd phi;		 // Volume porosity at position of cell walls
    vd rh;		 // Hydraulic radius 

    vd vx;		 // Vertex positions    
    vd vS;		 // Vertex cross-sectional area
    vd vSf;		 // Fluid-occupied cross-sectional area for vertex
    vd vSs;		 // Solid-occupied cross-sectional area for vertex
    vd vVf;		 // Fluid-occupied volume of cell
    vd vVs;		 // Solid-occupied volume of cell
    
    vd vphi;			// Volume porosity of vertex
    vd vrh;			// Hydraulic radius of vertex
  private:
    void Celldata();		// Compute all cell data
    bool prismatic=false;
  };                            /* class Geom */
}                               // namespace segment

#endif /* _GEOM_H_ */
