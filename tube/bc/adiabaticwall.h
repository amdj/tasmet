#pragma once
#ifndef _ADIABATICWALL_H_
#define _ADIABATICWALL_H__

#include "tubebcvertex.h"
#include "stateeq.h"


namespace tube{

  class RightAdiabaticWall:public TubeBcVertex // Right adiabatichermal wall boundary
  {
  public:
    virtual void show() const;
    virtual string getType() const {return string("RightAdiabaticWall");}
    virtual TubeBcVertex* copy() const {return new RightAdiabaticWall(*this);}
    virtual enum connectpos connectPos() const {return connectpos::right;}
    virtual void initTubeVertex(us i,const Tube& thisseg);
    virtual const variable::var& pR() const {return pr;};
    virtual void setpR(const variable::var& o){pr=o;}
  protected:
    variable::var pr;
    StateR sr;
    
  };
  class LeftAdiabaticWall:public TubeBcVertex // Right adiabatichermal wall boundary
  {
  public:
    virtual void show() const;
    virtual string getType() const {return string("LeftAdiabaticWall");}
    virtual TubeBcVertex* copy() const {return new LeftAdiabaticWall(*this);}
    virtual enum connectpos connectPos() const {return connectpos::left;}
    virtual void initTubeVertex(us i,const Tube& thisseg);
  };


} // namespace tube


#endif /* _ADIABATICWALL_H_ */




