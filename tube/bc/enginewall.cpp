#include "tube.h"
#include "tasystem.h"
#include "enginewall.h"

// file: bcvertex.cpp, created March 20th, 2014.
// Author: J.A. de Jong
#include "twimpedance.h"
#include "tubevertex.h"


namespace tube{

  vd EngineWallContinuity::error(const TubeVertex& v) const {
    TRACE(10,"EngineWallContinuity::errror()");
    vd error=Continuity::error(v);
    assert(v.gc!=NULL);
    assert(v.gc->getSys()!=NULL);
    error(0)=v.gc->getSys()->getCurrentMass()-v.gc->getMass();
    return error;
  }
  dmat EngineWallContinuity::drhoi(const TubeVertex& v) const {
    TRACE(10,"EngineWallContinuity::drhoi()");
    dmat drhoi=Continuity::drhoi(v);
    drhoi.row(0).zeros();
    drhoi(0,0)=v.lg.vVf*5e2;
    return drhoi;
  }
  dmat EngineWallContinuity::drhoim1(const TubeVertex& v) const {
    TRACE(10,"EngineWallContinuity::drhoim1()");
    dmat drhoim1=Continuity::drhoim1(v);
    drhoim1.row(0).zeros();
    return drhoim1;
  }
  dmat EngineWallContinuity::drhoip1(const TubeVertex& v) const {
    TRACE(10,"EngineWallContinuity::drhoip1()");
    dmat drhoip1=Continuity::drhoip1(v);
    drhoip1.row(0).zeros();
    return drhoip1;
  }
  
  void EngineWall::initTubeVertex(us i,const Tube& thisseg)
  {
    TRACE(8,"EngineWall::Init(), vertex "<< i <<".");
    TubeVertex::initTubeVertex(i,thisseg);
    ewc.init(thisseg);
    eqs[0]=&ewc;
  }

  void EngineWall::show() const {
    cout << "EngineWall boundary condition. Applies over-all mass conservation" << "\n";
    TubeVertex::show();
  }

} // namespace tube












