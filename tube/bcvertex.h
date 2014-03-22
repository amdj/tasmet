// file: tubebc.h, created March 20th, 2014
// Author: J.A. de Jong

// bcvertex.h: external boundary conditions for tubes. This file
// contains the implementation of typical external boundary conditions
// for tubes as a custom vertex. Examples are adiabatic walls, isothermal walls and an
// adiabatic open pressure boundary conditions.

#include "tube.h"
#include "vertex.h"

namespace tube{
  class RightVertex:public TubeVertex{};
  class RightAdiabaticWallBcVertex:public RightVertex{ // Adiabatic wall boundary condition applied to the rightmost node of a tube
  public:
    RightAdiabaticWallBcVertex();
    ~RightAdiabaticWallBcVertex();
  };

} // namespace tube

