#include "tubebccell.h"
#include "var.h"
#include "jacobian.h"
#include "weightfactors.h"


namespace tube{
  using tasystem::JacRow;
  using tasystem::JacCol;

  TubeBcCell::TubeBcCell(us i,const Tube& t):
    Cell(i,t)
  {}
   vd TubeBcCell::extrapolateQuant(physquant q) const {
    TRACE(10,"TubeBcCell::extrapolateQuant()");
    switch(q){
    case massFlow:
      return extrapolateMassFlow();
      break;
    case momentumFlow:
      return extrapolateMomentumFlow();
      break;
    case energyFlow:
      // return e.extrapolateEnergyFlow();
      break;
    case heatFlow:
      return e.extrapolateHeatFlow();
      break;
    case solidHeatFlow:
      return se.extrapolateHeatFlow();
    case rhoRT:
      return extrapolateRhoRT();
    default:
      break;
    }
    WARN("Fatal: No physquant selected");
    abort();
    return vd2();               // To avoid compiler warning
  }
  JacRow TubeBcCell::dExtrapolateQuant(physquant q) const{
    TRACE(10,"TubeBcCell::dExtrapolateQuant()");
    switch(q){
    case massFlow:
      return dExtrapolateMassFlow();
      break;
    case momentumFlow:
      return dExtrapolateMomentumFlow();
      break;
    case energyFlow:
      // return e.dExtrapolateEnergyFlow();      break;
    case heatFlow:
      return e.dExtrapolateHeatFlow();
      break;
    case solidHeatFlow:
      return se.dExtrapolateHeatFlow();
      break;
    case rhoRT:
      return dExtrapolateRhoRT();
    default:
      break;
    }
    WARN("Fatal: No physquant selected");
    abort();
    return JacRow(0);           // To avoid compiler warning
  }
}                // namespace tube

