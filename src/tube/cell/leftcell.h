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
    // virtual const tasystem::var& pL() const {return pL_;}

    const d& Sfbc() const {return Sfl;}
    // OVerloaded virtuals from BcCell
    void show(us detailnr=1) const;
    const tasystem::var& mbc() const {return ml_;}
    const tasystem::var& TL() const { return Tbc_;}
  };  

} // namespace tube

#endif /* _LEFTCELL_H_ */


