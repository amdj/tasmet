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
    variable::var Tbc_,mHbc_;
  public:
    BcCell(us i,const Tube& t);
    virtual ~BcCell(){}
    virtual Pos getPos() const=0;
    virtual void init(const Cell* left,const Cell* right);
    vd extrapolateQuant(Physquant) const;
    tasystem::JacRow dExtrapolateQuant(Physquant) const;
    // virtual d getValueBc(Varnr,us freqnr) const=0;
    Equation* Eq(EqType et) {return eqs.at(et);}

    // Fluid cross-sectional area at bc
    virtual  const d& Sfbc() const=0;

    // Return mass flow at the cell wall
    virtual const variable::var& mbc() const=0;
    // Return temperature at the cell wall
    virtual const variable::var& Tbc() const {return Tbc_;}
    // Return total enthalpy flow at the cell wall
    virtual const variable::var& mHbc() const {return mHbc_;}
  protected:
    
  };

} // namespace tube

#endif /* _BCCELL_H_ */
