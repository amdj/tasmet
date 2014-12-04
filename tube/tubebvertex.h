#pragma once
#ifndef _TUBEBVERTEX_H_
#define _TUBEBVERTEX_H_

#include "tubevertex.h"


namespace tube{
  // enum connectpos{ left,right};	// Where to connect the boundary condition.
  
  class Tube;

  class LeftTubeVertex:public TubeVertex{
    variable::var TL_,pL_,TsL_,UL_,rhoL_;
  public:
    LeftTubeVertex(us i,const Tube& t);
    virtual void init(const TubeVertex* left,const TubeVertex* right);
    virtual ~LeftTubeVertex(){}
    virtual const variable::var& UL() const {return UL_;}
    virtual const variable::var& pL() const {return pL_;}
    virtual const variable::var& TL() const {return TL_;}
    virtual const variable::var& TsL() const {return TsL_;}
    virtual void show() const;
};  
  class RightTubeVertex:public TubeVertex{
    variable::var TR_,pR_,TsR_,UR_,rhoR_;
  public:
    RightTubeVertex(us i,const Tube& t);
    virtual ~RightTubeVertex(){}
    virtual void init(const TubeVertex* left,const TubeVertex* right);
    virtual const variable::var& rhoR() const {return rhoR_;}
    virtual const variable::var& UR() const {return UR_;}
    virtual const variable::var& pR() const {return pR_;}
    virtual const variable::var& TR() const {return TR_;}
    virtual const variable::var& TsR() const {return TsR_;}
    virtual void show() const;
};  

} // namespace tube

#endif /* _TUBEBVERTEX_H_ */
