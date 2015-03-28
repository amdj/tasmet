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
  #ifdef SWIG
  %feature("notabstract") HopkinsLaminarDuct;
  #endif
  class HopkinsLaminarDuct:public LaminarDuct{
    HopkinsHeatSource hopkinsheat;
    d Tl,Tr;
    vd dTwdx;
    bool Tset=false;
  public:
    vd Tmirror;
    HopkinsLaminarDuct(const Geom& geom,d Tl);
    HopkinsLaminarDuct(const Geom& geom,d Tl,d Tr);
    HopkinsLaminarDuct(const HopkinsLaminarDuct& o);
    HopkinsLaminarDuct& operator=(const HopkinsLaminarDuct&)=delete;
    virtual segment::Seg* copy() const{return new HopkinsLaminarDuct(*this);}

    #ifndef SWIG
    virtual void init(const tasystem::TaSystem& gc);
    virtual const HeatSource& getHeatSource() const { return hopkinsheat;}
    #endif
  };


  
} /* namespace tube */

#endif	// HOPKINSLAMINARDUCT_H_






