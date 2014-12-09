#pragma once
#ifndef _TUBEBVERTEX_H_
#define _TUBEBVERTEX_H_

#include "tubevertex.h"
#include "varnr.h"

namespace segment{
  class Connector;
}
namespace tasystem{
  class JacRow;
}
namespace tube{

  enum pos{ left,right};	// Where to connect the boundary condition.
  class Tube;

  class TubeBcVertex{
    segment::Connector* con;
  public:
    virtual pos getPos() const=0;
    virtual ~TubeBcVertex(){}
    vd extrapolateQuant(physquant) const=0;
    JacRow dExtrapolateQuant(physquant) const=0;
  };
  class LeftTubeVertex:public TubeVertex,TubeBcVertex{

  public:
    
    LeftTubeVertex(us i,const Tube& t);
    virtual void init(const TubeVertex* left,const TubeVertex* right);
    virtual ~LeftTubeVertex(){}

    virtual pos getPos() const {return pos::left;}

    virtual const variable::var& UL() const;
    virtual const variable::var& pL() const;
    virtual const variable::var& TL() const;
    virtual const variable::var& TsL() const;
    virtual void show(us detailnr=1) const;
  
  };  

  class RightTubeVertex:public TubeVertex{
  public:
    RightTubeVertex(us i,const Tube& t);
    virtual ~RightTubeVertex(){}
    virtual void init(const TubeVertex* left,const TubeVertex* right);
    virtual pos getPos() const {return pos::right;}


    virtual const variable::var& rhoR() const;
    virtual const variable::var& UR() const;
    virtual const variable::var& pR() const;
    virtual const variable::var& TR() const;
    virtual const variable::var& TsR() const;
    virtual void show(us detailnr=1) const;

  };  

} // namespace tube

#endif /* _TUBEBVERTEX_H_ */
