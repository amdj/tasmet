// file: bcvertex.h, created March 20th, 2014
// Author: J.A. de Jong

// bcvertex.h: external boundary conditions for tubes. This file
// contains the implementation of typical external boundary conditions
// for tubes as a custom vertex. Examples are adiabatic walls, isothermal walls and an
// adiabatic open pressure boundary conditions.
#pragma once
#ifndef _TUBEBCVERTEX_H_
#define _TUBEBCVERTEX_H_


#include "tubevertex.h"


namespace tube{
  enum connectpos{ left,right};	// Where to connect the boundary condition.
  
  class TubeBcVertex:public TubeVertex
  {
  public:
    virtual string getType() const=0;
    virtual enum connectpos connectPos() const=0;
    virtual TubeBcVertex* copy() const=0; // Return a copy casted as class
    // TubeBcVertex. Used when init(gc) is run.
  };
  

} // namespace tube

#endif /* _TUBEBCVERTEX_H_ */
