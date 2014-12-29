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
  private:
    virtual vd extrapolateMassFlow() const=0;
    virtual tasystem::JacRow dExtrapolateMassFlow() const=0;
    virtual vd extrapolateDensity() const=0;
    virtual tasystem::JacRow dExtrapolateDensity() const=0;
    virtual vd extrapolateMomentumFlow() const=0;
    virtual tasystem::JacRow dExtrapolateMomentumFlow() const=0;
  };

} // namespace tube

#endif /* _TUBEBVERTEX_H_ */
