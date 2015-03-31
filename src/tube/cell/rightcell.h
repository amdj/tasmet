#pragma once
#ifndef _RIGHTCELL_H_
#define _RIGHTCELL_H_
#include "constants.h"
#include "bccell.h"
#include "state.h"

namespace tube{

  class RightCell:public BcCell{
    variable::var rhoUR_;
  public:
    RightCell(us i,const Tube& t);
    virtual ~RightCell(){}
    virtual void init(const Cell* left,const Cell* right);
    virtual Pos getPos() const {return Pos::right;}
    virtual const variable::var& rhoUR() const {return rhoUR_;}
    virtual void show(us detailnr=1) const;
    virtual d getValueBc(Varnr,us freqnr) const;
    virtual void setResVar(Varnr,const vd& res);
  private:
    virtual vd extrapolateMassFlow() const;
    virtual tasystem::JacRow dExtrapolateMassFlow() const;
    virtual vd extrapolateRhoRT() const;
    virtual tasystem::JacRow dExtrapolateRhoRT() const;
    virtual vd extrapolateMomentumFlow() const;
    virtual tasystem::JacRow dExtrapolateMomentumFlow() const;

  };  

} // namespace tube

#endif /* _RIGHTCELL_H_ */
