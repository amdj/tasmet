#pragma once
#ifndef _LEFTTUBEVERTEX_H_
#define _LEFTTUBEVERTEX_H_
#include "constants.h"
#include "tubebcvertex.h"

namespace tube{

  class LeftTubeVertex:public TubeBcVertex{
    variable::var UL_,TL_,TsL_,rhoL_;  
  public:
    LeftTubeVertex(us i,const Tube& t);
    virtual void init(const TubeVertex* left,const TubeVertex* right);
    virtual ~LeftTubeVertex(){}

    virtual Pos getPos() const {return Pos::left;}
    // For the pressure, we only do not assert that a vertex is
    // located on the left...
    virtual const variable::var& pL() const {return pL_;}
    virtual const variable::var& rhoL() const {return rhoL_;}
    virtual const variable::var& UL() const {return UL_;}
    virtual const variable::var& TL() const {return TL_;}
    virtual const variable::var& TsL() const {return TsL_;}
    virtual void show(us detailnr=1) const;

    virtual void setResVar(varnr,const vd& res);
    virtual d getValueBc(varnr,us freqnr) const;
  private:
    virtual vd extrapolateMassFlow() const;
    virtual tasystem::JacRow dExtrapolateMassFlow() const;
    virtual vd extrapolateRhoRT() const;
    virtual tasystem::JacRow dExtrapolateRhoRT() const;
    virtual vd extrapolateMomentumFlow() const;
    virtual tasystem::JacRow dExtrapolateMomentumFlow() const;
  };  

} // namespace tube

#endif /* _LEFTTUBEVERTEX_H_ */
