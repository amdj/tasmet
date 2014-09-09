#include "tube.h"
#include "adiabaticwall.h"

// file: bcvertex.cpp, created March 20th, 2014.
// Author: J.A. de Jong
#include "tubevertex.h"


namespace tube{
  void AdiabaticWall::show() const {
    cout << getType() << " boundary condition. Time-averaged part of prescribed temperature: " << Tbc << "\n";
    TubeVertex::show();
  }
  void AdiabaticWall::initTubeVertex(us i,const Tube& thisseg)
  {
    TRACE(8,"RightAdiabaticWall::Init(), vertex "<< i <<".");
    TubeVertex::initTubeVertex(i,thisseg);
    updateW(thisseg);
  }
  
  void RightAdiabaticWall::updateW(const SegBase& thisseg){
    TRACE(8,"RightAdiabaticWall::updateW()");
    
    eWc1=-w.vSfL/w.dxm;
    eWc2= w.vSf/w.dxm;
    eWc3= w.vSfR/lg.xr;
    eWc4=0;
    
  }
  vd RightAdiabaticWall::esource() const {
    TRACE(6,"RightAdiabaticWall::esource()");
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
  void LeftAdiabaticWall::updateW(const SegBase& thisseg){
    TRACE(8,"LeftAdiabaticWall::updateW()");
    
    eWc1=0;
    eWc2= w.vSfL/lg.xl;
    eWc3= w.vSf/w.dxp;
    eWc4=-w.vSfR/w.dxp;
    
  }
  vd LeftAdiabaticWall::esource() const {
    TRACE(6,"LeftAdiabaticWall::esource()");
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












