#include "righttubevertex.h"

namespace tube{
  using variable::var;
  using tasystem::JacRow;
  using tasystem::JacCol;

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

    // Initialize right wall variables
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

} // namespace tube
