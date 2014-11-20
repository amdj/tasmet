#pragma once

#ifndef _GEOMHELPERS_H_
#define _GEOMHELPERS_H_
#include "geom.h"
inline double min(double x,double y)  {
  return x<=y? x : y;
}
inline double max(double x,double y)  {
  return x<=y? y : x;
}


namespace geom{
  const int FIRST=0;
  const int LAST=1;
  // Smooth geometries of segments together. The percentage is a
  // length part of the smallest of the two given geometries.
  void smoothEnds(Geom& smooththisone,int smooththisonepos,
                  const Geom& tothisone,int tothisonepos,int perc);
  
  class VertPlates:public Geom{
    d S_=-1,phi_=1,rh_=-1;
  public:
    VertPlates(const Grid& g,d S,d phi,d y0,bool blapprox=false);
    VertPlates(const VertPlates& t):VertPlates(t.grid(),t.S_,t.phi_,t.rh_,t.isBlApprox()){} 
    virtual Geom* copy() const {return new VertPlates(*this);}
    virtual void show() const;
    string shape() const {return "vert";}
    virtual d S(us i) const {return S_;}		 // Cross sectional area as a function of x
    virtual d phi(us i) const {return phi_;}		 // Volume porosity at position of cell walls
    virtual d rh(us i) const {return rh_;}		 // Hydraulic radius
    virtual ~VertPlates(){}
  };


 

} // namespace segment

#endif /* _GEOMHELPERS_H_ */
