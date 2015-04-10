// isentropictube.h
//
// Author: J.A. de Jong 
//
// Description:
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef ISENTROPICTUBE_H
#define ISENTROPICTUBE_H

#include "tube.h"
#include "vtypes.h"
#include "math_common.h"

namespace tube{

  class IsentropicTube:public Tube
  {
    DragResistance nodrag;
    HeatSource noheat;
    IsentropicTube& operator=(const IsentropicTube&) =delete;    
    IsentropicTube(const IsentropicTube&)=delete;
    IsentropicTube(const IsentropicTube&,const tasystem::TaSystem&);
  public:
    IsentropicTube(const Geom& geom);
    virtual segment::Seg* copy(const tasystem::TaSystem&) const;
    virtual ~IsentropicTube();
    #ifndef SWIG
    vd dragCoefVec(us i) const;
    virtual const DragResistance& getDragResistance() const {return nodrag;}
    virtual const HeatSource& getHeatSource() const {return noheat;}
    void cleanup();
    #endif  // ifndef SWIG
  };
  
} /* namespace tube */


#endif // ISENTROPICTUBE_H
//////////////////////////////////////////////////////////////////////
