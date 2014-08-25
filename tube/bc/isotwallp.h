#pragma once
#ifndef _ISOTWALLP_H_
#define _ISOTWALLP_H_
#include "isotwall.h"
#include "momentumeq.h"
namespace tube{

  class IsoTWallPMomentum:public Momentum{
  public:
    virtual vd error(const TubeVertex&) const;
    virtual dmat dpi(const TubeVertex&) const;
    virtual dmat dpip1(const TubeVertex&) const;
    virtual dmat dpim1(const TubeVertex&) const;
    d mWp0im1,mWp0i,mWp0ip1;    
  };
  
  class LeftIsoTWallP:public LeftIsoTWall{
    d p0;
    IsoTWallPMomentum isotpmom;    
  public:
    LeftIsoTWallP(d Tbc,d p0=0):LeftIsoTWall(Tbc),p0(p0){}
    LeftIsoTWallP(const LeftIsoTWallP& o):LeftIsoTWallP(o.Tbc,o.p0){}
    LeftIsoTWallP& operator=(const LeftIsoTWallP& o){p0=o.p0; LeftIsoTWall::operator=(o);}
    virtual TubeBcVertex* copy() const {return new LeftIsoTWallP(*this);}
    virtual void initTubeVertex(us i,const Tube& t);
    virtual vd msource() const;
  protected:
    virtual void updateW(const SegBase&);
  };



} // namespace tube

#endif /* _ISOTWALLP_H_ */
