#include "tube.h"
#include "isotwall.h"

// file: bcvertex.cpp, created March 20th, 2014.
// Author: J.A. de Jong
#include "twimpedance.h"
#include "tubevertex.h"


namespace tube{
  void IsoTWall::show() const {
    cout << getType() << " boundary condition. Time-averaged part of prescribed temperature: " << Tbc << "\n";
    TubeVertex::show();
  }
  void IsoTWall::initTubeVertex(us i,const Tube& thisseg)
  {
    TRACE(8,"RightIsoTWall::Init(), vertex "<< i <<".");
    TubeVertex::initTubeVertex(i,thisseg);
    updateW(thisseg);
  }
  
  void RightIsoTWall::updateW(const SegBase& thisseg){
    TRACE(8,"RightIsoTWall::updateW()");
    
    eWc1=-w.vSfL/w.dxm;
    eWc2= w.vSf/w.dxm;
    eWc3= w.vSfR/lg.xr;
    eWc4=0;
    
  }
  vd RightIsoTWall::esource() const {
    TRACE(6,"RightIsoTWall::esource()");
    const dmat& fDFT=gc->fDFT;
    const dmat& iDFT=gc->iDFT;
    // Source term related to temperature boundary condition
    vd esource(gc->Ns,fillwith::zeros);
    TRACE(0,"Tbc:"<<Tbc);
    vd TRt=Tbc*vd(gc->Ns,fillwith::ones);
    vd T0=gc->T0*vd(gc->Ns,fillwith::ones);
    vd kappaR=static_cast<const Energy*>(eqs[2])->kappaR(*this);    
    // vd kappaR=gc->gas.kappa(T0);

    esource+=-1.0*w.vSfR*fDFT*(kappaR%TRt)/lg.xr;
    TRACE(3,"esource:"<<esource);
    return esource;  
  }
  void LeftIsoTWall::updateW(const SegBase& thisseg){
    TRACE(8,"LeftIsoTWall::updateW()");
    
    eWc1=0;
    eWc2= w.vSfL/lg.xl;
    eWc3= w.vSf/w.dxp;
    eWc4=-w.vSfR/w.dxp;
    
  }
  vd LeftIsoTWall::esource() const {
    TRACE(6,"LeftIsoTWall::esource()");
    const dmat& fDFT=gc->fDFT;
    const dmat& iDFT=gc->iDFT;
    // Source term related to temperature boundary condition
    vd esource(gc->Ns,fillwith::zeros);
    TRACE(0,"Tbc:"<<Tbc);
    vd TLt=Tbc*vd(gc->Ns,fillwith::ones);
    vd T0=gc->T0*vd(gc->Ns,fillwith::ones);
    vd kappaL=static_cast<const Energy*>(eqs[2])->kappaL(*this);    
    // vd kappaR=gc->gas.kappa(T0);
    esource+=-1.0*w.vSfL*fDFT*(kappaL%TLt)/lg.xl;
    TRACE(3,"esource:"<<esource);
    return esource;  
  }
} // namespace tube












