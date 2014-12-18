#include "lefttubevertex.h"

namespace tube{
  using variable::var;
  using tasystem::JacRow;
  using tasystem::JacCol;

  LeftTubeVertex::LeftTubeVertex(us i,const Tube& t):
    TubeBcVertex(i,t)
  {
    TRACE(15,"LeftTubeVertex::LeftTubeVertex()");
    rhoL_=var(gc);
    UL_=var(gc);
    TL_=var(gc);
    TsL_=var(gc);
  }
  void LeftTubeVertex::init(const TubeVertex* left,const TubeVertex* right){
    TRACE(10,"LeftTubeVertex::init()");
    TubeVertex::init(left,right);
    assert(!left);
    assert(right);


    vars.push_back(&rhoL_);
    vars.push_back(&UL_);
    vars.push_back(&TL_);
    vars.push_back(&TsL_);

    // Initialize left wall variables
    rhoL_=rho();
    TL_=T();
    TsL_=Ts();
  }
  void LeftTubeVertex::show(us detailnr) const{
    cout << "------------- LeftTubeVertex ----------\n";
    TubeVertex::show(detailnr);
  }
  void LeftTubeVertex::setResVar(varnr v,const vd& res){
    TRACE(15,"LeftTubeVertex::setResVar()");
    switch(v){
    case varnr::rhoL:
      rhoL_.set(res);
      break;
    case varnr::TL:
      TL_.set(res);
      break;
    case varnr::UL:
      UL_.set(res);
      break;
    case varnr::TsL:
      TsL_.set(res);
      break;
    case varnr::pL:
      pL_.set(res);
      break;
    default:
      TubeVertex::setResVar(v,res);
      break;
    }
  }

} // namespace tube
