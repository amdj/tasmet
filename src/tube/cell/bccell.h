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
  protected:
    tasystem::var Tbc_,mHbc_;
  public:
    BcCell(us i,const Tube& t);
    virtual ~BcCell(){}
    virtual Pos getPos() const=0;
    virtual void init(const Cell* left,const Cell* right);
    vd extrapolateQuant(Varnr) const;
    tasystem::JacRow dExtrapolateQuant(Varnr) const;
    // virtual d getValueBc(Varnr,us freqnr) const=0;

    // Fluid cross-sectional area at bc
    virtual  const d& Sfbc() const=0;

    // Return mass flow at the cell wall
    virtual const tasystem::var& mbc() const=0;
    // Return temperature at the cell wall
    virtual const tasystem::var& Tbc() const {return Tbc_;}
    // Return total enthalpy flow at the cell wall
    virtual const tasystem::var& mHbc() const {return mHbc_;}
  protected:
    
  };

} // namespace tube

#endif /* _BCCELL_H_ */
