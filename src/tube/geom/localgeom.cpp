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

    Sfl=geom.Sf(i);
    Sfr=geom.Sf(i+1);
    Ssl=geom.Ss(i);
    Ssr=geom.Ss(i+1);

    xr=geom.x(i+1);		// Position of right cell wall
    xl=geom.x(i);			// Position of left cell wall
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
    cout <<"Sfl   :" << Sfl<<"\n";
    cout <<"Sfr   :" << Sfr<<"\n";    
    cout <<"vVf   :" << vVf<<"\n";
    cout <<"vrh   :" << vrh<<"\n";
    cout <<"xl    :" << xl<<"\n";            
    cout <<"xr    :" << xr<<"\n";

  }
} // namespace tube


