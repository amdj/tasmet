/*
 * tube.h
 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#pragma once
#ifndef HOPKINSLAMINARDUCT_H_
#define HOPKINSLAMINARDUCT_H_
#include "laminarduct.h"
#include "hopkinsheat.h"
namespace tube{

  class HopkinsLaminarDuct:public LaminarDuct{
    HopkinsHeatSource hopkinsheat;
    HeatSource noheatatall;
  public:
    HopkinsLaminarDuct(const Geom& geom):LaminarDuct(geom),hopkinsheat(*this){    type="HopkinsLaminarDuct";}
    HopkinsLaminarDuct(const HopkinsLaminarDuct& o);
    HopkinsLaminarDuct& operator=(const HopkinsLaminarDuct&);
    virtual const HeatSource& getHeatSource() const { return hopkinsheat;}
    // virtual const HeatSource& getHeatSource() const { return noheatatall;}
  };
  HopkinsLaminarDuct HopkinsLaminarDuctTs(const Geom& geom,d Ts);



  
} /* namespace tube */

#endif	// HOPKINSLAMINARDUCT_H_






