#pragma once
#ifndef _ISOTWALL_H_
#define _ISOTWALL_H__

#include "tubebcvertex.h"



namespace tube{

  class IsoTWall:public TubeBcVertex //  isothermal wall boundary condition
  {
  public:
    IsoTWall(d Tbc): Tbc(Tbc){}
    IsoTWall(const IsoTWall& o):IsoTWall(o.Tbc){}
    IsoTWall& operator=(const IsoTWall& o){Tbc=o.Tbc;}
    virtual ~IsoTWall(){}
    virtual void show() const;
    virtual void initTubeVertex(us i,const Tube& thisseg);

  protected:
    virtual void updateW(const SegBase&)=0;
    d Tbc;			// Temperature at boundary.
  };

  class RightIsoTWall:public IsoTWall // Right isothermal wall boundary
  {
  public:
    RightIsoTWall(d Tbc):IsoTWall(Tbc){}
    RightIsoTWall(const RightIsoTWall& o):RightIsoTWall(o.Tbc){}
    RightIsoTWall& operator=(const RightIsoTWall& o){Tbc=o.Tbc;}
    virtual string getType() const {return string("RightIsoTWall");}
    virtual TubeBcVertex* copy() const {return new RightIsoTWall(*this);}
    virtual enum connectpos connectPos() const {return connectpos::right;}
    virtual vd esource() const final;		// Source term for constant
					// temperature
  protected:
    virtual void updateW(const SegBase&);    
  };
  class LeftIsoTWall:public IsoTWall // Left isothermal wall boundary
  {
  public:
    LeftIsoTWall(d Tbc):IsoTWall(Tbc){}
    LeftIsoTWall(const LeftIsoTWall& o):LeftIsoTWall(o.Tbc){}
    LeftIsoTWall& operator=(const LeftIsoTWall& o){Tbc=o.Tbc;}
    virtual string getType() const {return string("LeftIsoTWall");}
    virtual TubeBcVertex* copy() const {return new LeftIsoTWall(*this);}
    virtual enum connectpos connectPos() const {return connectpos::left;}
    virtual vd esource() const final;		// Source term for constant
					// temperature
  protected:
    virtual void updateW(const SegBase&);    
  };


} // namespace tube


#endif /* _ISOTWALL_H_ */




