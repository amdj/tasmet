#include "tube.h"
#include "adiabaticwall.h"

// file: bcvertex.cpp, created March 20th, 2014.
// Author: J.A. de Jong
#include "tubevertex.h"


namespace tube{
  void AdiabaticWall::show() const {
    cout << getType() << " boundary condition.\n";
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
  void LeftAdiabaticWall::updateW(const SegBase& thisseg){
    TRACE(8,"LeftAdiabaticWall::updateW()");
    
    eWc1=0;
    eWc2= w.vSfL/lg.xl;
    eWc3= w.vSf/w.dxp;
    eWc4=-w.vSfR/w.dxp;
    
  }

} // namespace tube












