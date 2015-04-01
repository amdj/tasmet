#pragma once
#ifndef _RIGHTCELL_H_
#define _RIGHTCELL_H_
#include "constants.h"
#include "bccell.h"
#include "state.h"

namespace tube{

  class RightCell:public BcCell{
    variable::var mR_;
  public:
    RightCell(us i,const Tube& t);
    virtual ~RightCell(){}
    virtual void init(const Cell* left,const Cell* right);
    virtual Pos getPos() const {return Pos::right;}
    virtual void show(us detailnr=1) const;

    const variable::var& mbc() const{return mR_;}
    virtual const variable::var& mR() const { return mR_;}
    
    // OVerloaded virtuals from BcCell
    vd extrapolateQuant(Physquant) const;
    tasystem::JacRow dExtrapolateQuant(Physquant) const;


  };  

} // namespace tube

#endif /* _RIGHTCELL_H_ */
