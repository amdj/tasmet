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
    HopkinsLaminarDuct(const Geom& geom):LaminarDuct(geom),hopkinsheat(*this){}
    HopkinsLaminarDuct(const HopkinsLaminarDuct& o);
    HopkinsLaminarDuct& operator=(const HopkinsLaminarDuct&);
    virtual const HeatSource& getHeatSource() const { return hopkinsheat;}
    virtual SegBase* copy() const{return new HopkinsLaminarDuct(*this);}
    virtual string getType() const {return string("HopkinsLaminarDuct");}
  };
  HopkinsLaminarDuct HopkinsLaminarDuctTs(const Geom& geom,d Ts);
  HopkinsLaminarDuct HopkinsLaminarDuctTs(const Geom& geom,d Tl,d Tr);


  
} /* namespace tube */

#endif	// HOPKINSLAMINARDUCT_H_






