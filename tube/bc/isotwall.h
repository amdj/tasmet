#pragma once
#ifndef _ISOTWALL_H_
#define _ISOTWALL_H__

#include "adiabaticwall.h"



namespace tube{

  class RightIsoTWall:public RightAdiabaticWall // Right isothermal wall boundary
  {
    d Tbc;
  public:
    RightIsoTWall(d Tbc):Tbc(Tbc){}
    RightIsoTWall(const RightIsoTWall& o):RightIsoTWall(o.Tbc){}
    virtual void initTubeVertex(us i,const Tube&);
    virtual string getType() const {return string("RightIsoTWall");}
    virtual TubeBcVertex* copy() const {return new RightIsoTWall(*this);}
    virtual void show() const;
    virtual enum connectpos connectPos() const {return connectpos::right;}
    virtual vd esource() const final;		// Source term for constant
    // temperature
  private:
    void updateW(const SegBase&);    
  };


} // namespace tube


#endif /* _ISOTWALL_H_ */




