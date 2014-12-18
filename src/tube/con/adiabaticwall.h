#pragma once
#ifndef _ADIABATICWALL_H_
#define _ADIABATICWALL_H__

#include "tubebc.h"


namespace tube{

 // Adiabatic wall boundary
  class AdiabaticWall:public TubeBc {

  public:
    AdiabaticWall(const TubeVertex& v): TubeBc(v){}
    virtual void show() const;
    virtual string getType() const {return string("AdiabaticWall");}
    virtual Connector* copy() const {return new RightAdiabaticWall(*this);}
    virtual enum connectpos connectPos() const {return position;}
    virtual void init();
  };

} // namespace tube


#endif /* _ADIABATICWALL_H_ */




