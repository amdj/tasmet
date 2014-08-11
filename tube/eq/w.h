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
    d wLl,wRr,wLr,wRl;		// Basis weight functions
    d wL0,wL1,wRNm1,wRNm2;    	// Special boundary weight functions
    d vxim1,vxi,vxip1;
    d dxm,dxp;
    d vSfR,vSfL;
    int UsignL;
    int UsignR;
    W();
    void show() const;
    void operator()(const tube::TubeVertex& v);
  };   	// Class W
}	// namespace W
#endif /* _W_H_ */
