#include "tubebvertex.h"
#include "var.h"
namespace tube{
  using variable::var;

  LeftTubeVertex::LeftTubeVertex(us i,const Tube& t):
    TubeVertex(i,t)
  {
    TRACE(15,"LeftTubeVertex::LeftTubeVertex()");
    UL_=var(gc);
    TL_=var(gc);
    pL_=var(gc);
    TsL_=var(gc);

  }
  vd LeftTubeBcVertex::extrapolateQuant(physquant q) const {
    TRACE(10,"LeftTubeVertex::extrapolateQuant()");
    const WeightFactors& w=v.weightFactors();

    switch(q){
    case massFlow:
      return c.extrapolateMassFlow();
      break;
    }
  }

  JacRow LeftTubeVertex::dExtrapolateQuant(physquant q) const{
    TRACE(10,"LeftTubeVertex::dExtrapolateQuant()");
    const WeightFactors& w=v.weightFactors();
    switch(q){
    case massFlow:
      return c.dExtrapolateMassFlow();
      break;
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
    TubeVertex(i,t)
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

  }
  void RightTubeVertex::show(us detailnr) const{
    cout << "------------- RightTubeVertex ---------\n";
    TubeVertex::show(detailnr);
  }
}                // namespace tube

