// file: bccell.h, created March 20th, 2014
// Author: J.A. de Jong

// bccell.h: external boundary conditions for tubes. This file
// contains the implementation of typical external boundary conditions
// for tubes as a custom cell. Examples are adiabatic walls, isothermal walls and an
// adiabatic open pressure boundary conditions.
#pragma once

#ifndef _IMPEDANCEBC_H_
#define _IMPEDANCEBC_H_

#include "constants.h"
#include "tubebc.h"

namespace tube{
  SPOILNAMESPACE

  class ImpedanceBc:public TubeBc {
  public:
    ImpedanceBc(us segnr,Pos pos,const tasystem::var& z,d T0=constants::T0);
  virtual ~ImpedanceBc();
};



} // namespace tube


#endif /* _IMPEDANCEBC_H_ */



