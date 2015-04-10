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

namespace tasystem{
  class TaSystem;
}
namespace tube{
  #ifdef SWIG
  %feature("notabstract") HopkinsLaminarDuct;
  #endif
  class HopkinsLaminarDuct:public LaminarDuct{
    HopkinsHeatSource hopkinsheat;
    d Tl,Tr;
    vd dTwdx;
    bool Tset=false;
    HopkinsLaminarDuct(const HopkinsLaminarDuct& o)=delete;
    HopkinsLaminarDuct& operator=(const HopkinsLaminarDuct&)=delete;
  public:
    vd Tmirror;
    HopkinsLaminarDuct(const HopkinsLaminarDuct& o,const tasystem::TaSystem&);
    HopkinsLaminarDuct(const Geom& geom,d Tl);
    HopkinsLaminarDuct(const Geom& geom,d Tl,d Tr);
    segment::Seg* copy(const tasystem::TaSystem&) const;
    #ifndef SWIG
    virtual const HeatSource& getHeatSource() const { return hopkinsheat;}
    #endif
  };


  
} /* namespace tube */

#endif	// HOPKINSLAMINARDUCT_H_






