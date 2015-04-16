#include "utils.h"
#include "vtypes.h"
#include "bccell.h"
#include "var.h"
#include "tubeequation.h"
#include "continuity.h"
#include "momentum.h"
#include "energy.h"
#include "jacrow.h"

#define iDFT (v.gc->iDFT)
#define fDFT (v.gc->fDFT)
#define DDTfd (v.gc->DDTfd)
#define Ns (v.gc->Ns())

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
    Tbc_.setadata(0,gc->T0());
    mHbc_=var(*gc);

    vars.push_back(&Tbc_);
    vars.push_back(&mHbc_);
  }
  vd BcCell::extrapolateQuant(Physquant p) const {
    TRACE(5,"LeftCell::extrapolateQuant()");
    switch(p){
    case Physquant::MassFlow:
      return Continuity::extrapolateMassFlow(*this);            
    case Physquant::MomentumFlow:
      return Momentum::extrapolateMomentumFlow(*this);      
      break;
    case Physquant::Pressure:
      return Momentum::extrapolatePressure(*this);      
      break;
    case Physquant::HeatFlow:
      return Energy::extrapolateHeatFlow(*this);
      break;
    case Physquant::EnthalpyFlow:
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
    case Physquant::MassFlow:
      return Continuity::dExtrapolateMassFlow(*this);            
    case Physquant::MomentumFlow:
      return Momentum::dExtrapolateMomentumFlow(*this);      
      break;
    case Physquant::Pressure:
      return Momentum::dExtrapolatePressure(*this);      
      break;
    case Physquant::HeatFlow:
      return Energy::dExtrapolateHeatFlow(*this);
      break;
    case Physquant::EnthalpyFlow:
      return Energy::dExtrapolateEnthalpyFlow(*this);
      break;
    default:
      WARN("This is not yet implemented!");
      assert(false);

    }
  }

}                // namespace tube

