#include "tube.h"
#include "adiabaticwall.h"

// file: bcvertex.cpp, created March 20th, 2014.
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
    eqs.push_back(std::unique_ptr<TubeEquation>(sr.copy()));
  }
  
 

} // namespace tube












