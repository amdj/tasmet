#include "tubebcvertex.h"
#include "var.h"
#include "jacobian.h"

namespace tube{
  using tasystem::JacRow;

  TubeBcVertex::TubeBcVertex(us i,const Tube& t):
    TubeVertex(i,t)
  {}
   vd TubeBcVertex::extrapolateQuant(physquant q) const {
    TRACE(10,"TubeBcVertex::extrapolateQuant()");
    switch(q){
    case massFlow:
      return c.extrapolateMassFlow();
      break;
    case momentumFlow:
      return m.extrapolateMomentumFlow();
      break;
    case energyFlow:
      // return e.extrapolateEnergyFlow();
      break;
    case heatFlow:
      return e.extrapolateHeatFlow();
      break;
    case solidHeatFlow:
      return se.extrapolateHeatFlow();
    default:
      break;
    }
    WARN("Fatal: No physquant selected");
    abort();
    return vd2();               // To avoid compiler warning
  }
  JacRow TubeBcVertex::dExtrapolateQuant(physquant q) const{
    TRACE(10,"TubeBcVertex::dExtrapolateQuant()");
    switch(q){
    case massFlow:
      return c.dExtrapolateMassFlow();
      break;
    case momentumFlow:
      return m.dExtrapolateMomentumFlow();
      break;
    case energyFlow:
      // return e.dExtrapolateEnergyFlow();      break;
    case heatFlow:
      return e.dExtrapolateHeatFlow();
      break;
    case solidHeatFlow:
      return se.dExtrapolateHeatFlow();
      break;
    default:
      break;
    }
    WARN("Fatal: No physquant selected");
    abort();
    return JacRow(0);           // To avoid compiler warning
  }
  
}                // namespace tube

