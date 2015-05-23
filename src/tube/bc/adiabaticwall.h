#pragma once
#ifndef _ADIABATICWALL_H_
#define _ADIABATICWALL_H__

#include "tubebccell.h"
#include "stateeq.h"


namespace tube{

  class RightAdiabaticWall:public BcCell // Right adiabatichermal wall boundary
  {
  public:
    virtual void show() const;
    virtual string getType() const {return string("RightAdiabaticWall");}
    virtual BcCell* copy() const {return new RightAdiabaticWall(*this);}
    virtual enum connectpos connectPos() const {return connectpos::right;}
    virtual void initCell(us i,const Tube& thisseg);
    virtual const tasystem::var& pR() const {return pr;};
    virtual void setpR(const tasystem::var& o){pr=o;}
  protected:
    tasystem::var pr;
    StateR sr;
    
  };
  class LeftAdiabaticWall:public BcCell // Right adiabatichermal wall boundary
  {
  public:
    virtual void show() const;
    virtual string getType() const {return string("LeftAdiabaticWall");}
    virtual BcCell* copy() const {return new LeftAdiabaticWall(*this);}
    virtual enum connectpos connectPos() const {return connectpos::left;}
    virtual void initCell(us i,const Tube& thisseg);
  };


} // namespace tube


#endif /* _ADIABATICWALL_H_ */




