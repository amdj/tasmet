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
#include <vtypes.h>
#include <math_common.h>
#include "continuityeq.h"
#include "momentumeq.h"
#include "isentropiceq.h"
#include "stateeq.h"
#include "solidenergyeq.h"


namespace tube{
  SPOILNAMESPACE

  using namespace segment;

  class IsentropicTube:public Tube
  {
  public:
    IsentropicTube(Geom geom);
    IsentropicTube(const IsentropicTube&);
    IsentropicTube& operator=(const IsentropicTube&);
    virtual void init(const Globalconf&);
    void cleanup();
    vector<const TubeEquation*> getEq() const;
    Continuity c;		// Continuity equation
    Momentum m;			// Momentum equation
    Isentropic is;		// Isentropic energy
    State s;			// State equation (ideal gas)
    Solidenergy se;		// Solid energy equation
    virtual ~IsentropicTube();
  };
  
} /* namespace tube */

#endif /* ISENTROPICTUBE_H_ */






