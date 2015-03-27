#pragma once
#ifndef _BOUNDARYLAYER_H_
#define _BOUNDARYLAYER_H_
#include <exception>
#include "vtypes.h"

namespace tube{
  #ifndef SWIG
  SPOILNAMESPACE;
  #endif

  class BoundaryLayer{
  protected:
    d dxb;                      // Smallest element
    d alpha;                    // Growth factor
    us n;                       // number of gridpoints in boundary layer
    BoundaryLayer(){}
  public:
    BoundaryLayer(d dxb,d L,d alpha);
    BoundaryLayer(d dxb,d L,us n); // Smallest, largest, number of gridpoints
    #ifndef SWIG
    vd getx() const;
    #endif
    virtual BoundaryLayer* copy() const {return new BoundaryLayer(*this);}
    virtual ~BoundaryLayer(){}
  };

  #ifndef SWIG
  class Grid;
  #endif

  // Special boundary layer such that grid spacing at the end matches
  // with grid spacing of Grid class given
  class AutoBoundaryLayer:public BoundaryLayer{
  public:
    AutoBoundaryLayer(d dxb,d alpha,const Grid& g) throw(std::exception);
    virtual BoundaryLayer* copy() const {return new AutoBoundaryLayer(*this);}
};


} // namespace tube


#endif /* _BOUNDARYLAYER_H_ */
