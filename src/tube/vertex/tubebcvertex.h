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

} // namespace tube

#endif /* _TUBEBVERTEX_H_ */
