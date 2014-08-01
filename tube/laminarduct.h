/*
 * tube.h
 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#pragma once
#ifndef LAMINARDUCT_H_
#define LAMINARDUCT_H__
#include "tube.h"
#include <vtypes.h>
#include <math_common.h>
#include "continuityeq.h"
#include "momentumeq.h"
#include "energyeq.h"
#include "stateeq.h"
#include "solidenergyeq.h"


namespace tube{
  SPOILNAMESPACE

  using namespace segment;

  class LaminarDuct:public Tube
  {
  public:
    LaminarDuct(Geom geom);
    LaminarDuct(const LaminarDuct&);
    LaminarDuct& operator=(const LaminarDuct&);
    virtual void init(const Globalconf&);
    void cleanup();
    vector<const TubeEquation*> getEq() const;
    Continuity c;		// Continuity equation
    Momentum m;			// Momentum equation
Energy e;			// Full energy equation
    State s;			// State equation (ideal gas)
    Solidenergy se;		// Solid energy equation
    virtual ~LaminarDuct();
  };
  
} /* namespace tube */

#endif	// LAMINARDUCT_H_






