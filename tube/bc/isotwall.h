#pragma once
#ifndef _ISOTWALL_H_
#define _ISOTWALL_H__

#include "tubebcvertex.h"



namespace tube{
  using segment::connectpos;

  class RightIsoTWall:public TubeBcVertex // Right isothermal wall boundary
  {
  public:
    RightIsoTWall(us segnr,d Tbc);
    RightIsoTWall(const RightIsoTWall& o);
    RightIsoTWall& operator=(const RightIsoTWall&);
    ~RightIsoTWall(){}
    virtual vd esource() const;		// Source term for constant temperature
    virtual void initTubeVertex(us i,const Tube& thisseg);
    virtual string gettype() const {return string("RightIsoTWall");}
    virtual enum connectpos connectPos() const {return connectpos::right;}
    d xhalf;			// Distance to right wall
  private:
    virtual void updateW(const SegBase&);
    d Tbc;
    
  };



} // namespace tube


#endif /* _ISOTWALL_H_ */



