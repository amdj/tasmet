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
#include "laminardrag.h"

namespace tube{
  SPOILNAMESPACE

  class LaminarDuct:public Tube
  {
    LaminarDragResistance laminardrag;
    HeatSource heat;
    LaminarDuct& operator=(const LaminarDuct&);
  public:
    LaminarDuct(const Geom& geom);
    LaminarDuct(const LaminarDuct&);
    virtual bool init(const tasystem::TaSystem&);
    virtual const DragResistance& getDragResistance() const {return laminardrag;}
    virtual const HeatSource& getHeatSource() const=0; // Yup,
    // abstract class. HopkinsLaminarTube implements  one version of
    // this heat source. Later on, AnneLaminarTube implements another
    // version of the heat source. This version is more applicable to
    // wide tubes for arbitrary cross-sectional geometries.
    void cleanup();
    virtual ~LaminarDuct();
    virtual vd dragCoefVec(us freqnr=1) const;
  };

  
} /* namespace tube */

#endif	// LAMINARDUCT_H_






