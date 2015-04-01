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
  public:
    BcCell(us i,const Tube& t);
    virtual ~BcCell(){}
    virtual Pos getPos() const=0;
    virtual vd extrapolateQuant(Physquant) const=0;
    virtual tasystem::JacRow dExtrapolateQuant(Physquant) const=0;
    // virtual d getValueBc(Varnr,us freqnr) const=0;
    Equation* Eq(us i) {return eqs.at(i);}
    virtual const variable::var& mbc() const=0; // Return mass flow at the
                                        // cell wall
  private:
  };

} // namespace tube

#endif /* _BCCELL_H_ */
