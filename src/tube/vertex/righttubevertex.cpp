#include "righttubevertex.h"

namespace tube{
  using variable::var;
  using tasystem::JacRow;
  using tasystem::JacCol;

  RightTubeVertex::RightTubeVertex(us i,const Tube& t):
    TubeBcVertex(i,t),
    sR(*this)
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
    sR.init();

    eqs.push_back(&sR);

    vars.push_back(&rhoR_);
    vars.push_back(&UR_);
    vars.push_back(&TR_);
    vars.push_back(&TsR_);
    vars.push_back(&pR_);    
    pR_=pL_;
    rhoR_=rho_;
    TR_=T_;
    TsR_=Ts_;
  }
  void RightTubeVertex::show(us detailnr) const{
    cout << "------------- RightTubeVertex ---------\n";
    TubeVertex::show(detailnr);
  }
  void RightTubeVertex::setResVar(varnr v,const vd& res){
    TRACE(15,"RightTubeVertex::setResVar()");
    switch(v){
    case varnr::rhoR:
      rhoR_.set(res);
      break;
    case varnr::TR:
      TR_.set(res);
      break;
    case varnr::UR:
      UR_.set(res);
      break;
    case varnr::TsR:
      TsR_.set(res);
      break;
    case varnr::pR:
      pR_.set(res);
      break;
    default:
      TubeVertex::setResVar(v,res);
      break;
    }
  }

} // namespace tube

