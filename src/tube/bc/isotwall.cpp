#include "tube.h"
#include "isotwall.h"

// file: bccell.cpp, created March 20th, 2014.
// Author: J.A. de Jong
#include "twimpedance.h"
#include "cell.h"


namespace tube{
  void RightIsoTWall::show() const {
    cout << getType() << " boundary condition. Time-averaged part of prescribed temperature: " << Tbc << "\n";
    Cell::show();
  }
  void RightIsoTWall::initCell(us i,const Tube& thisseg)
  {
    TRACE(8,"RightIsoTWall::initCell, cell "<< i <<".");
    RightAdiabaticWall::initCell(i,thisseg);
    RightIsoTWall::updateW(thisseg);
   }

   void RightIsoTWall::updateW(const SegBase& thisseg){
     TRACE(8,"RightIsoTWall::updateW()");

     e.Wc1=-lg.SfL/dxm;
     e.Wc2= lg.SfL/dxm;
     e.Wc3= lg.SfR/lg.xr;
     e.Wc4=0;

   }
   vd RightIsoTWall::esource() const {
     TRACE(6,"RightIsoTWall::esource()");
     const dmat& fDFT=gc->fDFT;
     const dmat& iDFT=gc->iDFT;
     // Source term related to temperature boundary condition
     vd esource(gc->Ns(),fillwith::zeros);
     TRACE(0,"Tbc:"<<Tbc);
     vd TRt=Tbc*vd(gc->Ns(),fillwith::ones);
     vd T0=gc->T0*vd(gc->Ns(),fillwith::ones);
     vd kappaR=e.kappaR(*this);    
     // vd kappaR=gc->gas.kappa(T0);
     esource+=-1.0*lg.SfR*fDFT*(kappaR%TRt)/lg.xr;
     TRACE(3,"esource:"<<esource);
     return esource;  
   }
  void LeftIsoTWall::show() const {
    cout << getType() << " boundary condition. Time-averaged part of prescribed temperature: " << Tbc << "\n";
    Cell::show();
  }
  void LeftIsoTWall::initCell(us i,const Tube& thisseg)
  {
    TRACE(8,"LeftIsoTWall::initCell, cell "<< i <<".");
    LeftAdiabaticWall::initCell(i,thisseg);
    LeftIsoTWall::updateW(thisseg);
   }

   void LeftIsoTWall::updateW(const SegBase& thisseg){
     TRACE(8,"LeftIsoTWall::updateW()");

     e.Wc1=0;
     e.Wc2= lg.SfL/lg.xl;
     e.Wc3= lg.SfR/dxp;
     e.Wc4=-lg.SfR/dxp;
   }
   vd LeftIsoTWall::esource() const {
     TRACE(6,"LeftIsoTWall::esource()");
     const dmat& fDFT=gc->fDFT;
     const dmat& iDFT=gc->iDFT;
     // Source term related to temperature boundary condition
     vd esource(gc->Ns(),fillwith::zeros);
     TRACE(0,"Tbc:"<<Tbc);
     vd TLt=Tbc*vd(gc->Ns(),fillwith::ones);
     vd T0=gc->T0*vd(gc->Ns(),fillwith::ones);
     vd kappaL=e.kappaL(*this);    
     // vd kappaR=gc->gas.kappa(T0);
     esource+=-1.0*lg.SfL*fDFT*(kappaL%TLt)/lg.xl;
     TRACE(3,"esource:"<<esource);
     return esource;  
   }

} // namespace tube












