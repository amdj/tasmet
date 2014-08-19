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
#include "vtypes.h"
#include "math_common.h"
#include "continuityeq.h"
#include "momentumeq.h"
#include "energyeq.h"
#include "stateeq.h"
#include "solidenergyeq.h"
#include "energyeq.h"
#include "laminardrag.h"

namespace tube{
  SPOILNAMESPACE


  using namespace segment;

  class LaminarDuct:public Tube
  {
  public:
    Continuity c;		// Continuity equation
    Momentum m;			// Momentum equation
    Energy e;			// Full energy equation
    State s;			// State equation (ideal gas)
    SolidTPrescribed se;		// Solid energy equation
  private:
    LaminarDragResistance laminardrag;
    HeatSource heat;

  public:
    LaminarDuct(const Geom& geom);
    LaminarDuct(const LaminarDuct&);
    LaminarDuct& operator=(const LaminarDuct&);
    virtual void init(const Globalconf&);
    virtual const DragResistance& getDragResistance() const {return laminardrag;}
    virtual const HeatSource& getHeatSource() const=0; // Yup,
    // abstract class. HopkinsLaminarTube implements  one version of
    // this heat source. Later on, AnneLaminarTube implements another
    // version of the heat source. This version is more applicable to
    // wide tubes for arbitrary cross-sectional geometries.
    void cleanup();
    EqVec getEq() const;
    virtual ~LaminarDuct();
  };

  
} /* namespace tube */

#endif	// LAMINARDUCT_H_






