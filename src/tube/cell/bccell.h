#pragma once
#ifndef _BCCELL_H_
#define _BCCELL_H_

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

  class BcCell:public Cell{
    segment::Connector* con=nullptr;
  public:
    BcCell(us i,const Tube& t);
    virtual ~BcCell(){}
    virtual Pos getPos() const=0;
    vd extrapolateQuant(physquant) const;
    tasystem::JacRow dExtrapolateQuant(physquant) const;
    virtual d getValueBc(Varnr,us freqnr) const=0;
  private:
    virtual vd extrapolateMassFlow() const=0;
    virtual tasystem::JacRow dExtrapolateMassFlow() const=0;
    virtual vd extrapolateRhoRT() const=0;
    virtual tasystem::JacRow dExtrapolateRhoRT() const=0;
    virtual vd extrapolateMomentumFlow() const=0;
    virtual tasystem::JacRow dExtrapolateMomentumFlow() const=0;
  };

} // namespace tube

#endif /* _BCCELL_H_ */
