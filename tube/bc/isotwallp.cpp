#include "isotwallp.h"

namespace tube{


  dmat IsoTWallPMomentum::dpi(const TubeVertex& v) const{
    TRACE(5,"IsoTWallPMomentum::dpi()");
    dmat dpi=Momentum::dpi(v);
    // TRACE(50,"dpi before zeroing"<< dpi);
    dpi.row(0).zeros();
    dpi(0,0)=mWp0i;
    // TRACE(50,"dpi after zeroing"<< dpi);
    return dpi;
  }
  dmat IsoTWallPMomentum::dpip1(const TubeVertex& v) const{
    TRACE(5,"IsoTWallPMomentum::dpip1()");
    dmat dpip1=Momentum::dpip1(v);
    // TRACE(50,"dpip1 before zeroing"<< dpip1);
    dpip1.row(0).zeros();
    dpip1(0,0)=mWp0ip1;
    // TRACE(50,"dpip1 after zeroing"<< dpip1);
    return dpip1;
  }
  dmat IsoTWallPMomentum::dpim1(const TubeVertex& v) const{
    TRACE(5,"IsoTWallPMomentum::dpim1()");
    dmat dpim1=Momentum::dpim1(v);
    dpim1.row(0).zeros();
    dpim1(0,0)=mWp0im1;
    return dpim1;
  }
  vd IsoTWallPMomentum::error(const TubeVertex& v) const{
    TRACE(5,"IsoTWallPMomentum::error()");
    vd error=Momentum::error(v);
    // Subtract old pressure contributions
    // and add new contributions
    error(0)+=(mWp0i-v.mWpi)*v.p(0);
    if(v.left!=NULL){
      // TRACE(30,"Adding left...");
      error(0)+=(mWp0im1-v.mWpim1)*v.left->p(0);
    }
    if (v.right!=NULL){
      // TRACE(30,"Adding right...");
      error(0)+=(mWp0ip1-v.mWpip1)*v.right->p(0);
    }
    // Source term is already added!
    return error;
  }
  void LeftIsoTWallP::initTubeVertex(us i,const Tube& thisseg){
    TRACE(50,"LeftIsoTWallP::initTubeVertex()");
    IsoTWall::initTubeVertex(i,thisseg);
    isotpmom.init(thisseg);
    eqs[1]=&isotpmom;
    
  }
  vd LeftIsoTWallP::msource() const{
    TRACE(5,"LeftIsoTWallPMomentum::msource()");
    vd msource(gc->Ns,fillwith::zeros);
    msource(0)=-1.0*w.vSfL*p0;
    return msource;
  }
  
  void LeftIsoTWallP::updateW(const SegBase& seg) {
    TRACE(10,"LeftIsoTWaLLP::updateW()");
    LeftIsoTWall::updateW(seg);
    isotpmom.mWp0im1=0;     
    isotpmom.mWp0i=w.wRl*w.vSfR+w.vSfL-w.vSfR;
    isotpmom.mWp0ip1=w.vSfR*w.wRr;
    // TRACE(20,"mWp0i:"<<isotpmom.mWp0i<<"\n");
    // TRACE(20,"mWp0ip1:"<<isotpmom.mWp0ip1<<"\n");  }

  }
} // namespace tube
