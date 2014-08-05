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
#include "energyeq_full.h"

namespace tube{
  SPOILNAMESPACE

  typedef vector<const TubeEquation*> EqVec;  
  using namespace segment;

  class LaminarDuct:public Tube
  {
  public:
    LaminarDuct(Geom geom);
    LaminarDuct(const LaminarDuct&);
    LaminarDuct& operator=(const LaminarDuct&);
    virtual void init(const Globalconf&);
    void cleanup();
    EqVec getEq() const;
    Continuity c;		// Continuity equation
    Momentum m;			// Momentum equation
    Energy e;			// Full energy equation
    State s;			// State equation (ideal gas)
    Solidenergy se;		// Solid energy equation
    virtual ~LaminarDuct();
  };

  class LaminarDuct_e: public LaminarDuct{
  public:
    LaminarDuct_e(Geom geom):LaminarDuct(geom) {  type="LaminarDuct_e";}
    LaminarDuct_e(const LaminarDuct_e& other): LaminarDuct_e(other.geom){}
    LaminarDuct_e& operator=(const LaminarDuct_e& other);
    virtual ~LaminarDuct_e() {}
    EqVec getEq() const;
    Energy_full e_full;  
  };  
  
} /* namespace tube */

#endif	// LAMINARDUCT_H_






