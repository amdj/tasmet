/*
 * tube.h
 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#pragma once
#ifndef ISENTROPICTUBE_H_
#define ISENTROPICTUBE_H__
#include "tube.h"
#include "vtypes.h"
#include "math_common.h"

namespace tube{
  SPOILNAMESPACE

  class IsentropicTube:public Tube
  {
    DragResistance nodrag;
    HeatSource noheat;
    IsentropicTube& operator=(const IsentropicTube&);    
  public:
    IsentropicTube(const Geom& geom);
    IsentropicTube(const IsentropicTube&);
    virtual const DragResistance& getDragResistance() const {return nodrag;}
    virtual const HeatSource& getHeatSource() const {return noheat;}
    virtual string getName() const {return string("IsentropicTube");}
    virtual segment::Seg* copy() const {TRACE(10,"IsentropicTube copy()");return new IsentropicTube(*this);}
    virtual void init(const tasystem::TaSystem&);
    void cleanup();
    virtual ~IsentropicTube();
  };
  
} /* namespace tube */

#endif /* ISENTROPICTUBE_H_ */





