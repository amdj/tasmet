#pragma once
#ifndef _LEFTCELL_H_
#define _LEFTCELL_H_
#include "constants.h"
#include "bccell.h"

namespace tube{

  class LeftCell:public BcCell{
  public:
    LeftCell(us i,const Tube& t);
    virtual void init(const Cell* left,const Cell* right);
    virtual ~LeftCell(){}

    virtual Pos getPos() const {return Pos::left;}
    // For the pressure, we only do not assert that a cell is
    // located on the left...
    // virtual const variable::var& pL() const {return pL_;}

    // OVerloaded virtuals from BcCell
    virtual void show(us detailnr=1) const;
    const variable::var& mbc() const {return mL_;}
    // const variable::var& mHbc() const {return mHL_;}
    virtual const variable::var& TL() const { return Tbc_;}

    vd extrapolateQuant(Physquant) const;
    tasystem::JacRow dExtrapolateQuant(Physquant) const;
  };  

} // namespace tube

#endif /* _LEFTCELL_H_ */


