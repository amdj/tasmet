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
  void LeftAdiabaticWall::initTubeVertex(us i,const Tube& thisseg)
  {
    TRACE(8,"RightAdiabaticWall::Init(), vertex "<< i <<".");
    TubeVertex::initTubeVertex(i,thisseg);
    updateW(thisseg);
  }
  void LeftAdiabaticWall::updateW(const SegBase& thisseg){
    TRACE(8,"LeftAdiabaticWall::updateW()");
  }


  // RightAdiabaticWall
  void RightAdiabaticWall::initTubeVertex(us i,const Tube& thisseg)
  {
    TRACE(8,"RightAdiabaticWall::Init(), vertex "<< i <<".");
    TubeVertex::initTubeVertex(i,thisseg);
    pr=var(gc);
    vars.push_back(&pr);
  }
  
 

} // namespace tube












