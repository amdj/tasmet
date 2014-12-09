#pragma once
#ifndef _TUBEBVERTEX_H_
#define _TUBEBVERTEX_H_

#include "tubevertex.h"
#include "var.h"
#include "varnr.h"
#include "pos.h"

namespace segment{
  class Connector;
}
namespace tasystem{
  class JacRow;
}
namespace tube{

  class Tube;

  class TubeBcVertex:public TubeVertex{
    segment::Connector* con;
  public:
    TubeBcVertex(us i,const Tube& t);
    virtual ~TubeBcVertex(){}
    virtual segment::pos getPos() const=0;
    vd extrapolateQuant(physquant) const;
    tasystem::JacRow dExtrapolateQuant(physquant) const;
  };
  class LeftTubeVertex:public TubeBcVertex{
    variable::var UL_,TL_,TsL_,rhoL_;  
  public:
    LeftTubeVertex(us i,const Tube& t);
    virtual void init(const TubeVertex* left,const TubeVertex* right);
    virtual ~LeftTubeVertex(){}

    virtual segment::pos getPos() const {return segment::pos::left;}

    virtual const variable::var& rhoL() const {return rhoL_;}
    virtual const variable::var& UL() const {return UL_;}
    virtual const variable::var& TL() const {return TL_;}
    virtual const variable::var& TsL() const {return TsL_;}
    virtual void show(us detailnr=1) const;
  };  
  class RightTubeVertex:public TubeBcVertex{
    variable::var rhoR_,TR_,TsR_,UR_,pR_;
  public:
    RightTubeVertex(us i,const Tube& t);
    virtual ~RightTubeVertex(){}
    virtual void init(const TubeVertex* left,const TubeVertex* right);
    virtual segment::pos getPos() const {return segment::pos::right;}


    virtual const variable::var& rhoR() const {return rhoR_;}
    virtual const variable::var& UR() const {return UR_;}
    virtual const variable::var& pR() const {return pR_;}
    virtual const variable::var& TR() const {return TR_;}
    virtual const variable::var& TsR() const {return TsR_;}
    virtual void show(us detailnr=1) const;

  };  

} // namespace tube

#endif /* _TUBEBVERTEX_H_ */
