#pragma once
#ifndef _TUBECONNECTOR_H_
#define _TUBECONNECTOR_H_

#include "tubevertex.h"

namespace segment{
  class Connector;
}

namespace tube{
  // enum connectpos{ left,right};	// Where to connect the boundary condition.
  
  class Tube;

  class RightTubeVertex:public TubeVertex{
  public:
    RightTubeVertex(us i,const Tube& t);
    virtual ~RightTubeVertex(){}
    virtual void init(const TubeVertex* left,const TubeVertex* right);
    virtual const variable::var& rhoR() const;
    virtual const variable::var& UR() const;
    virtual const variable::var& pR() const;
    virtual const variable::var& TR() const;
    virtual const variable::var& TsR() const;
    virtual void show(us detailnr=1) const;
};  

} // namespace tube

#endif /* _TUBECONNECTOR_H_ */
