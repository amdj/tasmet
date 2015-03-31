#pragma once
#ifndef _LEFTCELL_H_
#define _LEFTCELL_H_
#include "constants.h"
#include "bccell.h"

namespace tube{

  class LeftCell:public BcCell{
    variable::var UL_,TL_,TsL_,rhoL_;  
  public:
    LeftCell(us i,const Tube& t);
    virtual void init(const Cell* left,const Cell* right);
    virtual ~LeftCell(){}

    virtual Pos getPos() const {return Pos::left;}
    // For the pressure, we only do not assert that a cell is
    // located on the left...
    // virtual const variable::var& pL() const {return pL_;}
    virtual const variable::var& rhoL() const {return rhoL_;}
    virtual const variable::var& UL() const {return UL_;}
    virtual const variable::var& TL() const {return TL_;}
    virtual const variable::var& TsL() const {return TsL_;}
    virtual void show(us detailnr=1) const;

    virtual void setResVar(Varnr,const vd& res);
    virtual d getValueBc(Varnr,us freqnr) const;
  private:
    virtual vd extrapolateMassFlow() const;
    virtual tasystem::JacRow dExtrapolateMassFlow() const;
    virtual vd extrapolateRhoRT() const;
    virtual tasystem::JacRow dExtrapolateRhoRT() const;
    virtual vd extrapolateMomentumFlow() const;
    virtual tasystem::JacRow dExtrapolateMomentumFlow() const;
  };  

} // namespace tube

#endif /* _LEFTCELL_H_ */
