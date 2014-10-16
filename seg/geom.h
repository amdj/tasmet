#pragma once			// 
#ifndef _GEOM_H_
#define _GEOM_H_
#include "vtypes.h"
#include "localgeom.h"

#define MAXGP 500


namespace segment {
  SPOILNAMESPACE
  class Grid;
  void smoothEnds(Geom& first,int firstpos,Geom& second,int secondpos);
  class Geom{
    void createGeom(const vd& x,const vd& S,const vd& phi,const vd& rh,const string& cshape);
  public:
    Geom(){}
    Geom(const vd& x,const vd& S,const vd& phi,const vd& rh,const string& cshape);
    Geom(const Grid& g,const vd& S,const vd& phi,const vd& rh,const string& cshape);    
    Geom& operator=(const Geom& other);
    static Geom VertPlates(const Grid& g,d S,d phi,d y0);
    static Geom VertPlates(us gp,d L,d S,d phi,d y0);
    
    static Geom CylinderBlApprox(const Grid& g,d r);
    static Geom CylinderBlApprox(us gp,d L,d r);
    
    static Geom Cylinder(const Grid& g,d r);
    static Geom Cylinder(us gp,d L,d r);
    
    static Geom Cone(const Grid& g,d r1,d r2); // Return a cone
    static Geom Cone(us gp,d L,d r1,d r2); // Return a cone

    static Geom ConeBlApprox(const Grid& g,d r1,d r2); // Return a cone    blapprox
    static Geom ConeBlApprox(us gp,d L,d r1,d r2); // Return a cone    blapprox

    static Geom PrisVertStack(const Grid& g,d S,d phi,d rh); // Prismatic vertical plates stack
    static Geom PrisVertStack(us gp,d L,d S,d phi,d rh); // Prismatic vertical plates stack

    LocalGeom localGeom(us i) const;	// Get a local geometry for a certain vertex
    d getFluidVolume() const;
    bool prismatic=false;
    us gp;         // Number of cell walls
    us nCells;	 // Number of cells is gp-1
    d L;		 //Length of the Segment
    vd x;		 // Position of *CELL WALLS*
    vd S;		 // Cross sectional area as a function of x
    vd phi;		 // Volume porosity at position of cell walls
    string shape;	 // Shape keyword: currently available: 'circ','vert','blapprox'

    // These variables are all created with running the function Celldata()
    vd Ss;		 // Solid-occupied cross-sectional area
    vd Sf;		 // Fluid-occupied cross-sectional area of cell-wall

    vd rh;		 // Hydraulic radius 

    vd xv;		 // Vertex positions    
    vd vS;		 // Vertex cross-sectional area
    vd vSf;		 // Fluid-occupied cross-sectional area for vertex
    vd vSs;		 // Solid-occupied cross-sectional area for vertex
    vd vVf;		 // Fluid-occupied volume of cell
    vd vVs;		 // Solid-occupied volume of cell
    
    vd vphi;			// Volume porosity of vertex
    vd vrh;			// Hydraulic radius of vertex

    void show() const;
    void setPrismatic(bool isprismatic){prismatic=isprismatic;}

  };                            /* class Geom */
  


  
}                               // namespace segment

#endif /* _GEOM_H_ */
