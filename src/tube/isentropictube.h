// isentropictube.h
//
// Author: J.A. de Jong 
//
// Description:
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef ISENTROPICTUBE_H
#define ISENTROPICTUBE_H 1

#include "tube.h"
#include "vtypes.h"
#include "math_common.h"

namespace tube{

  class IsentropicTube:public Tube
  {
    DragResistance nodrag;
    HeatSource noheat;
    IsentropicTube& operator=(const IsentropicTube&) =delete;    
  public:
    IsentropicTube(const Geom& geom);
    IsentropicTube(const IsentropicTube&);
    virtual segment::Seg* copy() const {TRACE(10,"IsentropicTube copy()");return new IsentropicTube(*this);}
    virtual ~IsentropicTube();
    virtual string getType() const {return "IsentropicTube";}
    #ifndef SWIG
    vd dragCoefVec(us i) const;
    virtual const DragResistance& getDragResistance() const {return nodrag;}
    virtual const HeatSource& getHeatSource() const {return noheat;}
    virtual void init(const tasystem::TaSystem&);
    void cleanup();
    #endif
  };
  
} /* namespace tube */


#endif // ISENTROPICTUBE_H
//////////////////////////////////////////////////////////////////////
