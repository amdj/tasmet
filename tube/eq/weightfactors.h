#pragma once
#ifndef _WEIGHTFACTORS_H_
#define _WEIGHTFACTORS_H_
#include "localgeom.h"


namespace tube{
  class TubeVertex;
  class WeightFactors:public LocalGeom{
  public:
    d wLl=0,wRr=0,wLr=0,wRl=0;		// Basic weight functions
    d wL0=0,wL1=0,wRNm1=0,wRNm2=0;    	// Special boundary weight functions
    d vxm1=0,vx=0,vxp1=0;
    d dxm=0,dxp=0;
    d vSfR=0,vSfL=0;		// Cross sectional area at vx-position
    // of left and right neighbours

    WeightFactors(const TubeVertex& v);
    virtual void show() const;
  };

}      // namespace tube
#endif /* _WEIGHTFACTORS_H_ */

