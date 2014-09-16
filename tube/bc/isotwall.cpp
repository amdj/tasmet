#include "tube.h"
#include "isotwall.h"

// file: bcvertex.cpp, created March 20th, 2014.
// Author: J.A. de Jong
#include "twimpedance.h"
#include "tubevertex.h"


namespace tube{
  void RightIsoTWall::show() const {
    cout << getType() << " boundary condition. Time-averaged part of prescribed temperature: " << Tbc << "\n";
    TubeVertex::show();
  }
  void RightIsoTWall::initTubeVertex(us i,const Tube& thisseg)
  {
    TRACE(8,"RightIsoTWall::Init(), vertex "<< i <<".");
    RightAdiabaticWall::initTubeVertex(i,thisseg);
    RightIsoTWall::updateW(thisseg);
   }

   void RightIsoTWall::updateW(const SegBase& thisseg){
     TRACE(8,"RightIsoTWall::updateW()");

     eWc1=-lg.SfL/w.dxm;
     eWc2= lg.SfL/w.dxm;
     eWc3= lg.SfR/lg.xr;
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
     vd kappaR=static_cast<const Energy*>(eqs.at(2).get())->kappaR(*this);    
     // vd kappaR=gc->gas.kappa(T0);
     esource+=-1.0*lg.SfR*fDFT*(kappaR%TRt)/lg.xr;
     TRACE(3,"esource:"<<esource);
     return esource;  
   }

} // namespace tube












