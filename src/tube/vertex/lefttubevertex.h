#pragma once
#ifndef _LEFTTUBEVERTEX_H_
#define _LEFTTUBEVERTEX_H_
#include "pos.h"
#include "tubebcvertex.h"

namespace tube{

  class LeftTubeVertex:public TubeBcVertex{
    variable::var UL_,TL_,TsL_,rhoL_;  
  public:
    LeftTubeVertex(us i,const Tube& t);
    virtual void init(const TubeVertex* left,const TubeVertex* right);
    virtual ~LeftTubeVertex(){}

    virtual pos getPos() const {return segment::pos::left;}

    virtual const variable::var& rhoL() const {return rhoL_;}
    virtual const variable::var& UL() const {return UL_;}
    virtual const variable::var& TL() const {return TL_;}
    virtual const variable::var& TsL() const {return TsL_;}
    virtual void show(us detailnr=1) const;
  };  

} // namespace tube

#endif /* _LEFTTUBEVERTEX_H_ */
