// file: bcvertex.h, created March 20th, 2014
// Author: J.A. de Jong

// bcvertex.h: external boundary conditions for tubes. This file
// contains the implementation of typical external boundary conditions
// for tubes as a custom vertex. Examples are adiabatic walls, isothermal walls and an
// adiabatic open pressure boundary conditions.
#pragma once
#ifndef _BCVERTEX_H_
#define _BCVERTEX_H_


#include "tube.h"
#include "var.h"

namespace tube{
  class TubeVertex;
  class TubeBcVertex:public TubeVertex
  {
  public:
    TubeBcVertex(const Tube&t,us vertexnr);
    virtual ~TubeBcVertex(){}
  };
  



} // namespace tube

#endif /* _BCVERTEX_H_ */
