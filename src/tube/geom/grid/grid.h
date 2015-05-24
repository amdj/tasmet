#pragma once
#ifndef _GRID_H_
#define _GRID_H_

#include "vtypes.h"
#include <assert.h>
#include <exception>



namespace tube{
  #ifndef SWIG
  SPOILNAMESPACE
  #endif
  #ifdef SWIG

  #endif

  class BoundaryLayer;
  class Grid{
    BoundaryLayer *bL=nullptr,*bR=nullptr;
    us gp;                      // Only changeable by constructor
    d L;                        // Only changeable by constructor
    vd x;
  public:
    Grid(us gp,d L) throw(std::exception);
    // Copy constructor can just copy all
    Grid(const Grid& g);
    ~Grid();

    #ifndef SWIG
    Grid& operator=(const Grid& g);
    #endif
    // In these functions, we give the number of gridpoints (n) in the
    // boundary layer, plus the size of the layer (Lbl)

    // left boundary
    void setLeftBl(const BoundaryLayer& blleft); // Exponential growing
    // right boundary
    void setRightBl(const BoundaryLayer& blright); // Exponential growing

    bool isLeftBl() const { return bL?true:false;}
    bool isRightBl() const { return bR?true:false;}
    d getL() const {return L;}
    const vd& getx() const{return x;}
    us getgp() const {return gp;}
  private:
    void makex();
    void makeLeftBl();
    void makeRightBl();
   };

 } // namespace tube

#endif /* _GRID_H_ */

