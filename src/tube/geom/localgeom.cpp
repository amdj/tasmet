#include "localgeom.h"
#include "geom.h"
#include "assert.h"
#include "cell.h"
#include "tube.h"

namespace tube{

  LocalGeom Geom::localGeom(us i) const{
    assert(i<nCells());
    return LocalGeom(*this,i);
  }
  LocalGeom::LocalGeom(const Cell& v):
    LocalGeom(v.getTube().geom(),v.geti())
  {}
  LocalGeom::LocalGeom(const Geom& geom,us i)
  {
    TRACE(15,"LocalGeom::LocalGeom()");

    this->geom=&geom; 		// Save a pointer to the geometry instance
    assert(this->geom);
    this->i=i;
    vx=geom.vx(i);
    // initialize distances to next node to zero

    SfL=geom.Sf(i);
    SfR=geom.Sf(i+1);
    SsL=geom.Ss(i);
    SsR=geom.Ss(i+1);

    xR=geom.x(i+1);		// Position of right cell wall
    xL=geom.x(i);			// Position of left cell wall
    // Left and right cross-sectional area

    // Geometric parameters
    vSf=geom.vSf(i);
    vSs=geom.vSs(i);
    vVf=geom.vVf(i);
    vVs=geom.vVs(i);
    vrh=geom.vrh(i);
  }
  us LocalGeom::nCells() const{return geom->nCells();}	
  void LocalGeom::show() const {
    cout <<"Showing LocalGeom data..\n";
    cout <<"i     :" << i<<"\n";
    cout <<"vx    :" << vx<<"\n";
    cout <<"vSf   :" << vSf<<"\n";
    cout <<"SfL   :" << SfL<<"\n";
    cout <<"SfR   :" << SfR<<"\n";    
    cout <<"vVf   :" << vVf<<"\n";
    cout <<"vrh   :" << vrh<<"\n";
    cout <<"xL    :" << xL<<"\n";            
    cout <<"xR    :" << xR<<"\n";

  }
} // namespace tube


