#pragma once
#ifndef _RIGHTTUBEVERTEX_H_
#define _RIGHTTUBEVERTEX_H_
#include "pos.h"
#include "tubebcvertex.h"

namespace tube{

  class RightTubeVertex:public TubeBcVertex{
    variable::var rhoR_,TR_,TsR_,UR_,pR_;
  public:
    RightTubeVertex(us i,const Tube& t);
    virtual ~RightTubeVertex(){}
    virtual void init(const TubeVertex* left,const TubeVertex* right);
    virtual pos getPos() const {return pos::right;}

    virtual const variable::var& rhoR() const {return rhoR_;}
    virtual const variable::var& UR() const {return UR_;}
    virtual const variable::var& pR() const {return pR_;}
    virtual const variable::var& TR() const {return TR_;}
    virtual const variable::var& TsR() const {return TsR_;}
    virtual void show(us detailnr=1) const;

  };  

} // namespace tube

#endif /* _RIGHTTUBEVERTEX_H_ */
