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
  enum connectpos{ left,right};
  
  class TubeBcVertex:public TubeVertex
  {
    us segnumber;

  public:
    TubeBcVertex(us segnr);
    TubeBcVertex(const TubeBcVertex& o);
    TubeBcVertex& operator=(const TubeBcVertex& o);
    virtual ~TubeBcVertex(){}
    virtual string gettype() const=0;
    virtual enum connectpos connectPos() const=0;
    us segNumber() const {return segnumber;}
    void setSegNumber(us nr){segnumber=nr;}
    virtual TubeBcVertex* copy() const=0; // Return a copy casted as class
				    // Vertex. Used when init(gc) is run
				    // for a tube
  };
  

} // namespace tube

#endif /* _TUBEBCVERTEX_H_ */
