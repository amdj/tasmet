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
  #ifndef SWIG
  SPOILNAMESPACE
  #endif

  class LaminarDuct:public Tube
  {
    LaminarDragResistance laminardrag;
  protected:
    LaminarDuct(const Geom& geom);
    LaminarDuct(const LaminarDuct&)=delete;
    LaminarDuct(const LaminarDuct&,const tasystem::TaSystem&);
    LaminarDuct& operator=(const LaminarDuct&)=delete;
  public:
    #ifndef SWIG
    virtual ~LaminarDuct();
    virtual const DragResistance& getDragResistance() const {return laminardrag;}
    #endif
  };

  
} /* namespace tube */

#endif	// LAMINARDUCT_H_






