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

namespace tube{

  class HopkinsLaminarDuct:public LaminarDuct{
  public:
    HopkinsLaminarDuct(const Geom& geom):LaminarDuct(geom){}
    HopkinsLaminarDuct(const HopkinsLaminarDuct& o):HopkinsLaminarDuct(o.geom){}
    HopkinsLaminarDuct& operator=(const HopkinsLaminarDuct&);
  };

  HopkinsLaminarDuct ColdHopkinsLaminarDuct(const Geom& geom	\
					    ,const Globalconf& gc);



  
} /* namespace tube */

#endif	// HOPKINSLAMINARDUCT_H_






