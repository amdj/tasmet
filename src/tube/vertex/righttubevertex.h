#pragma once
#ifndef _RIGHTTUBEVERTEX_H_
#define _RIGHTTUBEVERTEX_H_
#include "constants.h"
#include "tubebcvertex.h"
#include "state.h"

namespace tube{

  class RightTubeVertex:public TubeBcVertex{
    variable::var rhoR_,TR_,TsR_,UR_,pR_;
    StateR sR;
  public:
    RightTubeVertex(us i,const Tube& t);
    virtual ~RightTubeVertex(){}
    virtual void init(const TubeVertex* left,const TubeVertex* right);
    virtual Pos getPos() const {return Pos::right;}

    virtual const variable::var& rhoR() const {return rhoR_;}
    virtual const variable::var& UR() const {return UR_;}
    virtual const variable::var& pR() const {return pR_;}
    virtual const variable::var& TR() const {return TR_;}
    virtual const variable::var& TsR() const {return TsR_;}

    virtual void show(us detailnr=1) const;
    virtual d getValueBc(varnr,us freqnr) const;
    virtual void setResVar(varnr,const vd& res);
  private:
    virtual vd extrapolateMassFlow() const;
    virtual tasystem::JacRow dExtrapolateMassFlow() const;
    virtual vd extrapolateRhoRT() const;
    virtual tasystem::JacRow dExtrapolateRhoRT() const;
    virtual vd extrapolateMomentumFlow() const;
    virtual tasystem::JacRow dExtrapolateMomentumFlow() const;

  };  

} // namespace tube

#endif /* _RIGHTTUBEVERTEX_H_ */
