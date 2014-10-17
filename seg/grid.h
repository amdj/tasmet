#pragma once
#ifndef _GRID_H_
#define _GRID_H_

#include "vtypes.h"
#include <assert.h>

namespace segment{
  SPOILNAMESPACE

  enum { left,right} LR;
  // class BoundaryLayer{


  // };
  
  class Grid{
    vd x;
    bool blleft=false;
    bool blright=false;
  public:
    Grid(us gp,d L);
    ~Grid(){}
    void setLeftBl(d dxmin,d percentage,us n); // Exponential growing left
                                     // boundary
    void setRightBl(d dxmin,d percentage,us n); // Exponential growing right boundary    
    const d& getL() const {return x(x.size()-1);}
    us getgp() const {return x.size();}
    const vd& getx() const {return x;}
   };

 }

#endif /* _GRID_H_ */

