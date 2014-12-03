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
namespace tasystem{
  class Globalconf;
}

namespace tube{
  enum connectpos{ left,right};	// Where to connect the boundary condition.
  
  class LeftTubeVertex:public TubeVertex{
  public:
    LeftTubeVertex(us i,const tasystem::Globalconf& g);
    virtual ~LeftTubeVertex();
};  
  class RightTubeVertex:public TubeVertex{
  public:
    RightTubeVertex(us i,const tasystem::Globalconf& g);
    virtual ~RightTubeVertex();
};  

} // namespace tube

#endif /* _TUBEBCVERTEX_H_ */
