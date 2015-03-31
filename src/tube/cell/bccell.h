#pragma once
#ifndef _TUBEBCELL_H_
#define _TUBEBCELL_H_

#include "cell.h"
#include "var.h"
#include "constants.h"

namespace segment{
  class Connector;
}
namespace tasystem{
  class JacRow;
}
namespace tube{

  class Tube;

  class TubeBcCell:public Cell{
    segment::Connector* con=NULL;
  public:
    TubeBcCell(us i,const Tube& t);
    virtual ~TubeBcCell(){}
    virtual Pos getPos() const=0;
    vd extrapolateQuant(physquant) const;
    tasystem::JacRow dExtrapolateQuant(physquant) const;
    virtual d getValueBc(varnr,us freqnr) const=0;
  private:
    virtual vd extrapolateMassFlow() const=0;
    virtual tasystem::JacRow dExtrapolateMassFlow() const=0;
    virtual vd extrapolateRhoRT() const=0;
    virtual tasystem::JacRow dExtrapolateRhoRT() const=0;
    virtual vd extrapolateMomentumFlow() const=0;
    virtual tasystem::JacRow dExtrapolateMomentumFlow() const=0;
  };

} // namespace tube

#endif /* _TUBEBCELL_H_ */
