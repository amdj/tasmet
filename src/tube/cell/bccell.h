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
    tasystem::var Tbc_,Tsbc_,mubc_,mHbc_,pbc_,rhobc_,ubc_;
  public:
    BcCell(us i,const Tube& t);
    virtual ~BcCell(){}
    friend class Tube;
    virtual Pos getPos() const=0;
    virtual void init(const Cell* left,const Cell* right);
    vd extrapolateQuant(Varnr) const;
    tasystem::JacRow dExtrapolateQuant(Varnr) const;
    // virtual d getValueBc(Varnr,us freqnr) const=0;

    // Fluid cross-sectional area at bc
    virtual  const d& Sfbc() const=0;

    // Return mass flow at the cell wall
    virtual const tasystem::var& mbc() const=0;
    // Return velocity at the cell wall
    // virtual const tasystem::var& ubc() const {return ubc_;}
    // // Return pressure at the cell wall
    const tasystem::var& pbc() const {return pbc_;}
    // Return density at the cell wall
    const tasystem::var& rhobc() const {return rhobc_;}
    // Return temperature at the cell wall
    const tasystem::var& Tbc() const {return Tbc_;}
    // Return  velocity at the cell wall
    const tasystem::var& Tsbc() const {return Tsbc_;}
    // Return  velocity at the cell wall
    const tasystem::var& ubc() const {return ubc_;}
    // Return total enthalpy flow at the cell wall
    const tasystem::var& mHbc() const {return mHbc_;}
  protected:
    
  };

} // namespace tube

#endif /* _BCCELL_H_ */
