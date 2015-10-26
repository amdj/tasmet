#pragma once

#ifndef _GEOMHELPERS_H_
#define _GEOMHELPERS_H_
#include "geom.h"
#include "constants.h"


namespace duct{

  #ifndef SWIG
  const int FIRST=0;
  const int LAST=1;
  // Smooth geometries of segments together. The percentage is a
  // length part of the smallest of the two given geometries.
  void smoothEnds(Geom& smooththisone,int smooththisonepos,
                  const Geom& tothisone,int tothisonepos,int perc);
  #endif

  class VertPlates:public Geom{
  protected:
    d S_=-1,phi_=1,rh_=-1;
  public:
    VertPlates(const vd& grid,d S,d phi,d y0,bool blapprox=false);
    VertPlates(const VertPlates& t):
      Geom(t),
      S_(t.S_),
      phi_(t.phi_),
      rh_(t.rh_)
    {}

    virtual Geom* copy() const {return new VertPlates(*this);}
    virtual void show() const;
    string shape() const {return "vert";}
    virtual d S(us i) const {return S_;}		 // Cross sectional area as a function of x
    virtual d phi(us i) const {return phi_;}		 // Volume porosity at position of cell walls
    virtual d rh(us i) const {return rh_;}		 // Hydraulic radius
    virtual ~VertPlates(){}
  };

} // namespace duct

#endif /* _GEOMHELPERS_H_ */



