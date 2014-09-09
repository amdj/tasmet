#pragma once
#ifndef _ADIABATICWALL_H_
#define _ADIABATICWALL_H__

#include "tubebcvertex.h"



namespace tube{

  class AdiabaticWall:public TubeBcVertex //  adiabatichermal wall boundary condition
  {
  public:
    virtual void show() const;
    virtual void initTubeVertex(us i,const Tube& thisseg);
  protected:
    virtual void updateW(const SegBase&)=0;
  };

  class RightAdiabaticWall:public AdiabaticWall // Right adiabatichermal wall boundary
  {
  public:
    virtual string getType() const {return string("RightAdiabaticWall");}
    virtual TubeBcVertex* copy() const {return new RightAdiabaticWall(*this);}
    virtual enum connectpos connectPos() const {return connectpos::right;}
    virtual vd esource() const final;		// Source term for constant
					// temperature
    virtual void initTubeVertex(us i,const Tube& thisseg);
  protected:
    void updateW(const SegBase&);    
  };
  class LeftAdiabaticWall:public AdiabaticWall // Left adiabatichermal wall boundary
  {
  public:
    virtual string getType() const {return string("LeftAdiabaticWall");}
    virtual TubeBcVertex* copy() const {return new LeftAdiabaticWall(*this);}
    virtual enum connectpos connectPos() const {return connectpos::left;}
  protected:
    void updateW(const SegBase&);    
  };


} // namespace tube


#endif /* _ADIABATICWALL_H_ */




