#pragma once
#ifndef _GRID_H_
#define _GRID_H_

#include "vtypes.h"
#include <assert.h>

namespace geom{
  SPOILNAMESPACE


  class BoundaryLayer;
  class Grid{
    BoundaryLayer *bL=NULL,*bR=NULL;
    us gp;                      // Only changeable by constructor
    d L;                        // Only changeable by constructor
    vd x;
    void makex();
    void makeLeftBl();
    void makeRightBl();
  public:
    Grid(us gp,d L);
    Grid(d L,d dx):Grid((unsigned int) int(ceil(L/dx))+1,L){} // Counts gp based on ceil(L/dx)
    // Copy constructor can just copy all
    Grid(const Grid& g);
    Grid& operator=(const Grid& g);
    ~Grid();
    // In these functions, we give the number of gridpoints (n) in the
    // boundary layer, plus the size of the layer (Lbl)
    void setLeftBl(const BoundaryLayer& blleft); // Exponential growing
    // left boundary
    void setRightBl(const BoundaryLayer& blright); // Exponential growing
    // right boundary

    bool isLeftBl() const { return bL?true:false;}
    bool isRightBl() const { return bR?true:false;}
    d getL() const {return L;}
    const vd& getx() const{return x;}
    us getgp() const {return gp;}
   };

 } // namespace geom

#endif /* _GRID_H_ */

