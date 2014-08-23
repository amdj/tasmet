#include "isotwallp.h"

namespace tube{


  dmat IsoTWallPMomentum::dpi(const TubeVertex& v) const{
    dmat dpi=Momentum::dpi(v);
    // Zero out first row
    dpi.row(0).zeros();
    // Add term
    // dpi(0,0)=v.
    return dpi;
  }
  vd IsoTWallPMomentum::error(const TubeVertex& v) const{
    vd error=Momentum::error(v);
    // Subtract old pressure contributions
    error(0)-=v.mWpi*v.p(0);
    error(0)-=v.mWpip1*v.right->p(0);

    // Source term is already added!
    WARN("Not yet done here!");
    exit(1);
    // Add new contributions
  }
  vd LeftIsoTWallP::msource() const{
    TRACE(5,"LeftIsoTWallPMomentum::msource()");
    vd msource(gc->Ns,fillwith::zeros);
    msource(0)=-1.0*lg.SfL*p0;
    return msource;
  }
  
  void LeftIsoTWallP::updateW(const SegBase& seg) {
    TRACE(6,"LeftIsoTWaLLP::updateW()");
    LeftIsoTWall::updateW(seg);
    mWp0im1=0;     
    mWp0i=w.wRl*w.vSfR+w.vSfL-w.vSfR;
    mWp0ip1=w.vSfR*w.wRr;
    
  }


}
