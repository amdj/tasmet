#pragma once
#ifndef _TUBEBVERTEX_H_
#define _TUBEBVERTEX_H_

#include "tubevertex.h"

namespace tube{
  // enum connectpos{ left,right};	// Where to connect the boundary condition.
  
  class LeftTubeVertex:public TubeVertex{
  public:
    LeftTubeVertex(us i,const tasystem::Globalconf& g);
    virtual ~LeftTubeVertex();
};  
  class RightTubeVertex:public TubeVertex{
  public:
    RightTubeVertex(us i,const tasystem::Globalconf& g);
    virtual ~RightTubeVertex();
};  

} // namespace tube

#endif /* _TUBEBVERTEX_H_ */
