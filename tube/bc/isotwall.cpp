#include "tube.h"
#include "isotwall.h"

// file: bcvertex.cpp, created March 20th, 2014.
// Author: J.A. de Jong
#include "twimpedance.h"
#include "tubevertex.h"


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

  void RightIsoTWall::initTubeVertex(us i,const Tube& thisseg)
  {
    TRACE(8,"RightIsoTWall::Init(), vertex "<< i <<".");
    TubeVertex::initTubeVertex(i,thisseg);
    updateW(thisseg);
  }
  
  void RightIsoTWall::updateW(const SegBase& thisseg){
    TRACE(8,"RightIsoTWall::updateW()");

    xhalf=lg.xR-lg.vxi;
    TRACE(30,"xhalf:"<<xhalf);

    // eWjim1=0;
    // eWji=0;
    // eWjip1=0;
    
    eWc1=-lg.SfL/lg.dxm;
    eWc2= lg.SfL/lg.dxm;
    eWc3= lg.SfR/lg.xr;
    eWc4=0;
    
  }

  vd RightIsoTWall::esource() const {
    // Source term related to temperature boundary condition
    TRACE(6,"RightIsoTWall::esource()");
    const dmat& fDFT=gc->fDFT;
    const dmat& iDFT=gc->iDFT;
    vd esource(gc->Ns,fillwith::zeros);
    TRACE(0,"Tbc:"<<Tbc);
    vd TRt=Tbc*vd(gc->Ns,fillwith::ones);
    vd T0=gc->T0*vd(gc->Ns,fillwith::ones);
    vd kappaR=static_cast<const Energy*>(eq[2])->kappaR(*this);    
    // vd kappaR=gc->gas.kappa(T0);
    esource+=-1.0*lg.SfR*fDFT*(kappaR%TRt)/lg.xr;
    return esource;  
  }

} // namespace tube












