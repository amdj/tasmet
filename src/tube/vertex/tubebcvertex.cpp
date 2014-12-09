#include "tubebcvertex.h"
#include "var.h"
#include "jacobian.h"

namespace tube{
  using variable::var;
  using tasystem::JacRow;
  using tasystem::JacCol;

  TubeBcVertex::TubeBcVertex(us i,const Tube& t):
    TubeVertex(i,t)
  {}
  LeftTubeVertex::LeftTubeVertex(us i,const Tube& t):
    TubeBcVertex(i,t)
  {
    TRACE(15,"LeftTubeVertex::LeftTubeVertex()");
    rhoL_=var(gc);
    UL_=var(gc);
    TL_=var(gc);
    TsL_=var(gc);
  }
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
      return e.extrapolateEnergyFlow();
      break;
    case heatFlow:
      // return e.extrapolateHeatFlow();
      break;
    case solidHeatFlow:
      // return s.extrapolateSolidHeatFlow();
      break;
    default:
      WARN("Fatal: No physquant selected");
      abort();
      return vd2();
    }
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
      // return e.dExtrapolateEnergyFlow();
      break;
    case heatFlow:
      // return e.dExtrapolateHeatFlow();
      break;
    case solidHeatFlow:
      // return s.dExtrapolateSolidHeatFlow();
      break;
    default:
      WARN("Fatal: No physquant selected");
      abort();
      return JacRow(0);         
    }
  }
  void LeftTubeVertex::init(const TubeVertex* left,const TubeVertex* right){
    TubeVertex::init(left,right);
    assert(!left);
    assert(right);
    // Initialize left wall variables

  }
  void LeftTubeVertex::show(us detailnr) const{
    cout << "------------- LeftTubeVertex ----------\n";
    TubeVertex::show(detailnr);
  }

  RightTubeVertex::RightTubeVertex(us i,const Tube& t):
    TubeBcVertex(i,t)
  {
    // Initialize right wall variables
    rhoR_=var(gc);
    UR_=var(gc);
    TR_=var(gc);
    pR_=var(gc);
    TsR_=var(gc);

  }
  void RightTubeVertex::init(const TubeVertex* left,const TubeVertex* right){
    TRACE(15,"RightTubeVertex::init()");
    assert(!right);             // Otherwise, this is not the
                                // rightmost!
    assert(left);
    TubeVertex::init(left,right);
    vars.push_back(&rhoR_);
    vars.push_back(&UR_);
    vars.push_back(&TR_);
    vars.push_back(&TsR_);
    vars.push_back(&pR_);    
  }
  void RightTubeVertex::show(us detailnr) const{
    cout << "------------- RightTubeVertex ---------\n";
    TubeVertex::show(detailnr);
  }
}                // namespace tube

