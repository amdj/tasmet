#pragma once
#ifndef _ENGINEWALL_H_
#define _ENGINEWALL_H__

#include "tubebcvertex.h"
#include "continuityeq.h"
// Special boundary condition which computes the time-avg density
// error based on conservation of total mass in the system


// Total mass in system is given:
// Totalmass: x kg
// Compute current mass in segment with:
// curmas= getCurMass() (integrate density over total volume of
// system)
// curmass= sum rho_i * Vf_i
// the error is: curmas - Totalmass
// The derivative of the error to the current density is: Vfi



namespace tube{

  class EngineWallContinuity:public Continuity
  {
  public:
    virtual vd error(const TubeVertex&) const;
    virtual dmat drhoi(const TubeVertex&) const;
    virtual dmat drhoim1(const TubeVertex&) const;
    virtual dmat drhoip1(const TubeVertex&) const;
  };
  
  class EngineWall:public TubeBcVertex // Wall which preserves global mass
  {
    EngineWallContinuity ewc;

  public:
    virtual void show() const;
    virtual void initTubeVertex(us i,const Tube& thisseg);

  };

  class LeftEngineWall:public EngineWall{
  public:
    virtual enum connectpos connectPos() const {return connectpos::left;}
    virtual string getType() const {return string("LeftEngineWall");}
    virtual TubeBcVertex* copy() const {return new LeftEngineWall(*this);}
    
  };


} // namespace tube


#endif /* _ENGINEWALL_H_ */



