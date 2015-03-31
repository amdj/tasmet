#include "tube.h"
#include "adiabaticwall.h"

// Author: J.A. de Jong
#include "cell.h"


namespace tube{
  void RightAdiabaticWall::show() const {
    cout << getType() << " boundary condition.\n";
    Cell::show();
  }
  void RightAdiabaticWall::initCell(us i,const Tube& thisseg)
  {
    TRACE(8,"RightAdiabaticWall::Init(), cell "<< i <<".");
    Cell::initCell(i,thisseg);
    pr=var(gc);
    vars.push_back(&pr);
    eqs.push_back(&sr);
  }
  void LeftAdiabaticWall::initCell(us i,const Tube& thisseg)
  {
    TRACE(8,"LeftAdiabaticWall::initCell(), cell "<< i <<".");
    Cell::initCell(i,thisseg);
  }
  void LeftAdiabaticWall::show() const {
    cout << getType() << " boundary condition.\n";
    Cell::show();
  }

} // namespace tube


