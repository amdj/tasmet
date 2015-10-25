#pragma once
#ifndef _GRID_H_
#define _GRID_H_

#include "vtypes.h"
#include <assert.h>
#include <exception>



namespace duct{
  #ifndef SWIG
  SPOILNAMESPACE
  #endif
  #ifdef SWIG
  // %catc
  #endif

  class Grid{
  public:
    Grid& operator=(const Grid& g)=delete;
    virtual vd getx() const=0;
    virtual ~Grid(){}
  };

  class LinearGrid:public Grid{
    us gp;
    d L;			// Length of the Duct
  public:
    LinearGrid(us gp,d L);
    LinearGrid(const LinearGrid& g):LinearGrid(g.gp,g.L){}
    vd getx() const {return linspace(0,L,gp);}
  };

  // Boundary layer grid.
  // L: length of grid
  // dxb: boundary layer grid spacing (minimal grid spacing)
  // xb: boundary layer tickness
  // dxmid: spacing in the middle part
  

  class BlGrid:public Grid{
    d L,dxb,dxmid;
  public:
    BlGrid(d L,d dxb,d dxmid);
    vd getx() const;
  };

} // namespace duct

#endif /* _GRID_H_ */

