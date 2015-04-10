#include "bccell.h"
#include "var.h"
#include "continuity.h"
#include "momentum.h"
#include "energy.h"
#include "jacrow.h"

namespace tube{
  using tasystem::JacRow;
  using tasystem::JacCol;
  using variable::var;
  BcCell::BcCell(us i,const Tube& t):
    Cell(i,t)
  {}
  void BcCell::init(const Cell* left,const Cell* right){
    TRACE(10,"BcCell::init(const Cell* left,const Cell* right)");
    Cell::init(left,right);
    Tbc_=var(*gc);
    mHbc_=var(*gc);
    Tbc_.setadata(0,gc->T0());

    vars.push_back(&Tbc_);
    vars.push_back(&mHbc_);
  }
  vd BcCell::extrapolateQuant(Physquant p) const {
    TRACE(5,"LeftCell::extrapolateQuant()");
    switch(p){
    case Physquant::massFlow:
      return Continuity::extrapolateMassFlow(*this);            
    case Physquant::momentumFlow:
      return Momentum::extrapolateMomentumFlow(*this);      
      break;
    case Physquant::heatFlow:
      return Energy::extrapolateHeatFlow(*this);
      break;
    case Physquant::enthalpyFlow:
      return Energy::extrapolateEnthalpyFlow(*this);
      break;
    default:
      WARN("This is not yet implemented!");
      assert(false);
    }

  }
  JacRow BcCell::dExtrapolateQuant(Physquant p) const {
    TRACE(5,"LeftCell::dExtrapolateQuant()");
    switch(p){
    case Physquant::massFlow:
      return Continuity::dExtrapolateMassFlow(*this);            
    case Physquant::momentumFlow:
      return Momentum::dExtrapolateMomentumFlow(*this);      
      break;
    case Physquant::heatFlow:
      return Energy::dExtrapolateHeatFlow(*this);
      break;
    case Physquant::enthalpyFlow:
      return Energy::dExtrapolateEnthalpyFlow(*this);
      break;
    default:
      WARN("This is not yet implemented!");
      assert(false);

    }
  }

}                // namespace tube

