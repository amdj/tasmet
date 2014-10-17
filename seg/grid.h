#pragma once
#ifndef _GRID_H_
#define _GRID_H_

#include "vtypes.h"
#include <assert.h>

namespace segment{
  SPOILNAMESPACE

  class Grid{
    bool blleft=false;
    bool blright=false;
    d percL,percR;
    d dxbL,dxbR;
    us nL,nR;
    us gp;
    d L;
    void makeLeftBl(vd& x) const;
    void makeRightBl(vd& x) const;
  public:
    Grid(us gp,d L);
    ~Grid(){}
    void setLeftBl(d dxmin,d percentage,us n); // Exponential growing left
                                     // boundary
    void setRightBl(d dxmin,d percentage,us n); // Exponential growing
                                                // right boundary
    const d& getL() const {return L;}
    vd getx() const;
   };

 }

#endif /* _GRID_H_ */

