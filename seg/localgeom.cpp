#include "localgeom.h"
#include "geom.h"
#include "assert.h"

namespace segment{

  LocalGeom Geom::localGeom(us i) const{
    assert(i<nCells);
    return LocalGeom(*this,i);
  }

  LocalGeom::LocalGeom(const Geom& geom,us i)
  {
    this->geom=&geom; 		// Save a pointer to the geometry instance
    const vd& vx=geom.vx;
    this->i=i;
    nCells=geom.nCells;
    vxi=vx(i);
    vxip1=0;			// To be filled below
    vxim1=0;			// To be filled below
    // initialize distances to next node to zero
    dxm=dxp=0;
    if(i>0){
      vxim1=vx(i-1);
      dxm=vxi-vxim1;
    }
    // ****************************** Initalization of vxipm and dxpm
    if(i<nCells-1){
      vxip1=vx(i+1);
      dxp=vxip1-vxi;
    }

    xR=geom.x(i+1);		// Position of right cell wall
    xL=geom.x(i);			// Position of left cell wall
    // Left and right cross-sectional area
    SfL=geom.Sf(i);
    SfR=geom.Sf(i+1);
    // Geometric parameters
    vSf=geom.vSf(i);
    vSs=geom.vSs(i);
    vVf=geom.vVf(i);
    vVs=geom.vVs(i);
    vrh=geom.vrh(i);
    xr=xR-vxi;
    xl=vxi-xL;
    assert(xl>0); assert(xr>0);
  }
  void LocalGeom::show(){

    cout<< "LocalGeom show() needs to be created.\n";
  }
} // namespace segment


