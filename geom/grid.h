#pragma once
#ifndef _GRID_H_
#define _GRID_H_

#include "vtypes.h"
#include <assert.h>

namespace geom{
  SPOILNAMESPACE

  class Grid{
    bool blleft=false;
    bool blright=false;
    d percL,percR;
    d dxbL,dxbR;
    us nL,nR;
    us gp;                      // Only changeable by constructor
    d L;                        // Only changeable by constructor
    vd x;
    void makeLeftBl();
    void makeRightBl();
    void makex();
  public:
    Grid(us gp,d L);
    // Copy constructor can just copy all
    ~Grid(){}
    void setLeftBl(d dxmin,d percentage,us n); // Exponential growing
    // left boundary

    bool isLeftBl() const {return blleft;}
    bool isRightBl() const { return blright;}
    void setRightBl(d dxmin,d percentage,us n); // Exponential growing
    // right boundary
    d getL() const {return L;}
    const vd& getx() const{return x;}
   };

 } // namespace geom

#endif /* _GRID_H_ */

