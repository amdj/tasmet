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
  using tasystem::var;

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
  vd BcCell::extrapolateQuant(Varnr v) const {
    TRACE(5,"LeftCell::extrapolateQuant()");
    switch(v){
    case Varnr::m:
      return Continuity::extrapolateMassFlow(*this);            
    case Varnr::mu:
      return Momentum::extrapolateMomentumFlow(*this);      
      break;
    case Varnr::p:
      return Momentum::extrapolatePressure(*this);      
      break;
    case Varnr::Q:
      return Energy::extrapolateHeatFlow(*this);
      break;
    case Varnr::mH:
      return Energy::extrapolateEnthalpyFlow(*this);
      break;
    default:
      WARN("This is not yet implemented!");
      assert(false);
    }
  }
  JacRow BcCell::dExtrapolateQuant(Varnr p) const {
    TRACE(5,"LeftCell::dExtrapolateQuant()");
    switch(p){
    case Varnr::m:
      return Continuity::dExtrapolateMassFlow(*this);            
    case Varnr::mu:
      return Momentum::dExtrapolateMomentumFlow(*this);      
      break;
    case Varnr::p:
      return Momentum::dExtrapolatePressure(*this);      
      break;
    case Varnr::Q:
      return Energy::dExtrapolateHeatFlow(*this);
      break;
    case Varnr::mH:
      return Energy::dExtrapolateEnthalpyFlow(*this);
      break;
    default:
      WARN("This is not yet implemented!");
      assert(false);

    }
  }

}                // namespace tube

