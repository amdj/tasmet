#include "isotwall.h"

// file: bcvertex.cpp, created March 20th, 2014.
// Author: J.A. de Jong
#include "twimpedance.h"
#include "tubevertex.h"
#include "momscale.h"


namespace tube{
  RightIsoTWall::RightIsoTWall(us segnr,d Tbc):TubeBcVertex(segnr),Tbc(Tbc){
    TRACE(8,"RightIsoTWall constructor");
    // Change continuity equation for open boundary
  }
  RightIsoTWall::RightIsoTWall(const RightIsoTWall& other):RightIsoTWall(other.segNumber(),other.Tbc)
  {
    TRACE(8,"RightIsoTWall copy cc.");
  }
  RightIsoTWall& RightIsoTWall::operator=(const RightIsoTWall& o){
    TRACE(8,"RightIsoTWall copy assignment operator");
    setSegNumber(o.segNumber());
    Tbc=o.Tbc;
    return *this;
  }

  void RightIsoTWall::Init(us i,const SegBase& thisseg)
  {
    TRACE(8,"RightIsoTWall::Init(), vertex "<< i <<".");
    TubeVertex::Init(i,thisseg);
    updateW(thisseg);
  }
  
  void RightIsoTWall::updateW(const SegBase& thisseg){
    TRACE(8,"RightIsoTWall::updateW()");

    xhalf=xR-vxi;
    TRACE(30,"xhalf:"<<xhalf);
    e.Wc1=-SfL/dxm;
    e.Wc2= SfL/dxm;
    e.Wc3= SfR/xhalf;
    e.Wc4=0;
    
  }

  vd RightIsoTWall::esource() const {
    // Source term related to temperature boundary condition
    TRACE(6,"RightIsoTWall::esource()");
    const dmat& fDFT=gc->fDFT;
    const dmat& iDFT=gc->iDFT;
    vd esource(gc->Ns,fillwith::zeros);
    TRACE(0,"Tbc:"<<Tbc);
    vd TRt=Tbc*vd(gc->Ns,fillwith::ones);
    vd kappaR=gc->gas.kappa(TRt);
    esource+=-1.0*SfR*fDFT*(kappaR%TRt)/xhalf;
    return esource;  
  }

} // namespace tube












