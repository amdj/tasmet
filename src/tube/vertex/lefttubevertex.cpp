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
    TubeVertex::init(left,right);
    assert(!left);
    assert(right);

    // Initialize left wall variables
    vars.push_back(&rhoL_);
    vars.push_back(&UL_);
    vars.push_back(&TL_);
    vars.push_back(&TsL_);
  }
  void LeftTubeVertex::show(us detailnr) const{
    cout << "------------- LeftTubeVertex ----------\n";
    TubeVertex::show(detailnr);
  }

} // namespace tube
