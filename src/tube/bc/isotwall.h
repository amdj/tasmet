#pragma once
#ifndef _ISOTWALL_H_
#define _ISOTWALL_H__

#include "adiabaticwall.h"



namespace tube{

  class RightIsoTWall:public RightAdiabaticWall // Right isothermal wall boundary
  {
    d Tbc;
  public:
    virtual void show() const;
    RightIsoTWall(d Tbc):Tbc(Tbc){}
    RightIsoTWall(const RightIsoTWall& o):RightIsoTWall(o.Tbc){}
    virtual void initCell(us i,const Tube&);
    virtual string getType() const {return string("RightIsoTWall");}
    virtual BcCell* copy() const {return new RightIsoTWall(*this);}
    virtual vd esource() const final;		// Source term for constant
    // temperature
  private:
    void updateW(const SegBase&);    
  };
  class LeftIsoTWall:public LeftAdiabaticWall // Left isothermal wall boundary
  {
    d Tbc;
  public:
    virtual void show() const;
    LeftIsoTWall(d Tbc):Tbc(Tbc){}
    LeftIsoTWall(const LeftIsoTWall& o):LeftIsoTWall(o.Tbc){}
    virtual void initCell(us i,const Tube&);
    virtual string getType() const {return string("LeftIsoTWall");}
    virtual BcCell* copy() const {return new LeftIsoTWall(*this);}
    virtual vd esource() const final;		// Source term for constant
    // temperature
  private:
    void updateW(const SegBase&);    
  };


} // namespace tube


#endif /* _ISOTWALL_H_ */




