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
    DragResistance nodrag;
    HeatSource noheat;
    
  public:
    IsentropicTube(Geom geom);
    IsentropicTube(const IsentropicTube&);
    IsentropicTube& operator=(const IsentropicTube&);
    virtual const DragResistance& getDragResistance() const {return nodrag;}
    virtual const HeatSource& getHeatSource() const {return noheat;}
    virtual string getName() const {return string("Isentropic tube");}
    virtual SegBase* copy() const {TRACE(10,"IsentropicTube copy()");return new IsentropicTube(*this);}
    virtual void init(const Globalconf&);
    void cleanup();
    vector<const TubeEquation*> getEqs() const;
    Continuity c;		// Continuity equation
    Momentum m;			// Momentum equation
    Isentropic is;		// Isentropic energy
    State s;			// State equation (ideal gas)
    SolidTPrescribed se;		// Solid energy equation
    virtual ~IsentropicTube();
  };
  
} /* namespace tube */

#endif /* ISENTROPICTUBE_H_ */






