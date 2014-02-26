/*
 * tube.h
 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */
#pragma once

#ifndef TUBE_H_
#define TUBE_H_
#include "globalconf.h"
#include "vvar.h"
#include "common/vtypes.h"

namespace tube {

  class Geom{
  public:
    us gp;		 /* Numberof gridpoints */
    vd x;			/* Grid points vector */
    vd S;			/* Cross sectional area as a function of x */
  };				/* class Geom */





  class tube {
  protected:
    globalconf::Globalconf& gc;
  };












} /* namespace tube */



#endif /* TUBE_H_ */
