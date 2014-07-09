// file: bcvertex.h, created March 20th, 2014
// Author: J.A. de Jong

// bcvertex.h: external boundary conditions for tubes. This file
// contains the implementation of typical external boundary conditions
// for tubes as a custom vertex. Examples are adiabatic walls, isothermal walls and an
// adiabatic open pressure boundary conditions.
#pragma once
#ifndef _BCVERTEX_H_
#define _BCVERTEX_H_

#include <vtypes.h>

namespace segment{

  enum connectpos{ left,right};  
  
  
  class BcVertex
  {

  public:
    BcVertex(us segnr):segnumber(segnr){TRACE(100,"Node not yet initialized!");}

    virtual ~BcVertex(){}
    // virtual Vertex* copy(const SegBase&)=0; // Copy the boundary condition vertex
    virtual string gettype() const=0;
    virtual enum connectpos connectPos() const=0;
    us segNumber() const {return segnumber;}
  private:
    us segnumber;
    
  };
  



} // namespace tube

#endif /* _BCVERTEX_H_ */
