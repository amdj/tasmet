#pragma once
#ifndef _W_H_
#define _W_H_
#include "geom.h"
#include "segbase.h"

namespace tube{  class TubeVertex;}
namespace W{
  SPOILNAMESPACE
  using segment::LocalGeom;

  class W{
  public:
    d wLl=0,wRr=0,wLr=0,wRl=0;		// Basis weight functions
    d wL0=0,wL1=0,wRNm1=0,wRNm2=0;    	// Special boundary weight functions
    d xvim1=0,xvi=0,xvip1=0;
    d dxm=0,dxp=0;
    d vSfR=0,vSfL=0;		// Cross sectional area at x-position
				// of vertex.
    int UsignL=1;
    int UsignR=1;
    W();
    void show() const;
    void operator()(const tube::TubeVertex& v);
  };   	// Class W
}	// namespace W
#endif /* _W_H_ */
