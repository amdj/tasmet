#pragma once
#ifndef _WEIGHTFACTORS_H_
#define _WEIGHTFACTORS_H_
#include "localgeom.h"


namespace tube{
  class Cell;
  class WeightFactors:public LocalGeom{
  public:
    d wLl=0,wRr=0,wLr=0,wRl=0;		// Basic weight functions
    d wL0=0,wL1=0,wRNm1=0,wRNm2=0;    	// Special boundary weight functions
    d vxm1=0,vxp1=0;               // Position of vx of the nodes
                                        // left and right of the current
    // vx(i-1) and vx(i+1) respectively
    d vSfR=0,vSfL=0;		// Cross sectional area at vx-position
    // of left and right neighbours

    WeightFactors(const Cell& v);
    virtual void show() const;
  };

}      // namespace tube
#endif /* _WEIGHTFACTORS_H_ */

