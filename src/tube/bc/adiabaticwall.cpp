#include "tube.h"
#include "adiabaticwall.h"

// Author: J.A. de Jong
#include "tubevertex.h"


namespace tube{
  void RightAdiabaticWall::show() const {
    cout << getType() << " boundary condition.\n";
    TubeVertex::show();
  }
  void RightAdiabaticWall::initTubeVertex(us i,const Tube& thisseg)
  {
    TRACE(8,"RightAdiabaticWall::Init(), vertex "<< i <<".");
    TubeVertex::initTubeVertex(i,thisseg);
    pr=var(gc);
    vars.push_back(&pr);
    eqs.push_back(&sr);
  }
  void LeftAdiabaticWall::initTubeVertex(us i,const Tube& thisseg)
  {
    TRACE(8,"LeftAdiabaticWall::initTubeVertex(), vertex "<< i <<".");
    TubeVertex::initTubeVertex(i,thisseg);
  }
  void LeftAdiabaticWall::show() const {
    cout << getType() << " boundary condition.\n";
    TubeVertex::show();
  }

} // namespace tube

