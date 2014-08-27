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
    d Tl,Tr;
    bool Tset=false;
  public:
    HopkinsLaminarDuct(const Geom& geom,d Tl);
    HopkinsLaminarDuct(const Geom& geom,d Tl,d Tr);
    HopkinsLaminarDuct(const HopkinsLaminarDuct& o);
    HopkinsLaminarDuct& operator=(const HopkinsLaminarDuct&);
    virtual void init(const Globalconf& gc);
    virtual const HeatSource& getHeatSource() const { return hopkinsheat;}
    virtual SegBase* copy() const{return new HopkinsLaminarDuct(*this);}
    virtual string getName() const {return string("Hopkins laminar duct");}
  };


  
} /* namespace tube */

#endif	// HOPKINSLAMINARDUCT_H_






