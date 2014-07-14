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
#include "bcvertex.h"

namespace tube{

  class TubeBcVertex:public TubeVertex,public BcVertex
  {
  public:
    TubeBcVertex(us segnr): BcVertex(segnr){
      TRACE(8,"TubeBcVertex(segnr)");
    }
    TubeBcVertex(TubeBcVertex& o):TubeBcVertex(o.segNumber()) {TRACE(8,"TubeBcVertex copy cc.");}
    TubeBcVertex& operator=(const TubeBcVertex& o)
    {
      TRACE(8,"TubeBcVertex::operator=()");
      setSegNumber(o.segNumber());
      return *this;      
    }
    virtual ~TubeBcVertex(){}
    // virtual Vertex* copy(const SegBase&)=0; // Copy the boundary condition vertex
  };
  



} // namespace tube

#endif /* _TUBEBCVERTEX_H_ */
