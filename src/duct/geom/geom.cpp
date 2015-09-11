#include "geom.h"
#include "grid.h"
#include <assert.h>

namespace duct{

  
  Geom::Geom(const Grid& grid,bool blapprox,bool prismatic):
    grid_(grid),
    blapprox(blapprox),
    prismatic(prismatic)
  {}

  // void Geom::show() const{
  //   cout << "No effort done to show geom\n";
  // }
  d Geom::vx(us i) const{
    assert(i<nCells());
    return 0.5*(x(i+1)+x(i));
  }
  d Geom::vphi(us i) const{
    assert(i<nCells());
    return 0.5*(phi(i+1)+phi(i));
  }
  d Geom::vrh(us i) const{
    assert(i<nCells());
    return 0.5*(rh(i+1)+rh(i));
  }

  d Geom::vS(us i) const {
    assert(i<nCells());
    return 0.5*(S(i+1)+S(i));
  }
  d Geom::vSs(us i) const {
    assert(i<nCells());
    return (1-vphi(i))*vS(i);
  }
  d Geom::vSf(us i) const {
    assert(i<nCells());
    return vphi(i)*vS(i);
  }
  vd Geom::vSf_vec() const{
    return vS_vec()%vphi_vec();
  }
  vd Geom::vS_vec() const{
    vd vS(nCells());
    for(us i=0;i<nCells();i++)
      vS(i)=this->vS(i);
    return vS;
  }
  vd Geom::vrh_vec() const{
    vd vrh(nCells());
    for(us i=0;i<nCells();i++)
      vrh(i)=this->vrh(i);
    return vrh;
  }
  vd Geom::vphi_vec() const{
    vd vphi(nCells());
    for(us i=0;i<nCells();i++)
      vphi(i)=this->vphi(i);
    return vphi;
  }
  vd Geom::vx_vec() const{
    TRACE(10,"Geom::vx_vec()");
    VARTRACE(5,nCells());
    vd vx(nCells());
    for(us i=0;i<nCells();i++)
      vx(i)=this->vx(i);
    return vx;
  }

  d Geom::vVf(us i) const {
    assert(i<nCells());
    return vSf(i)*(x(i+1)-x(i));
  }
  d Geom::vVs(us i) const {
    assert(i<nCells());
    return vSs(i)*(x(i+1)-x(i));
  }

  d Geom::getFluidVolume() const {
    us nCells=this->nCells();
    d Vf=0;
    for(us i=0;i<nCells;i++)
      Vf+=vVf(i);
    return Vf;
  }
  d Geom::getSolidVolume() const {
    us nCells=this->nCells();
    d Vs=0;
    for(us i=0;i<nCells;i++)
      Vs+=vVs(i);
    return Vs;
  }


} // namespace duct




