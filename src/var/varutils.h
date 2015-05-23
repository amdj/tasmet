// varutils.h
//
// Author: J.A. de Jong 
//
// Description:
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef VARUTILS_H
#define VARUTILS_H
#include "globalconf.h"
#include "var.h"

namespace tasystem {
  
  var coldTemp(const var& pres); // Returns a variable with constant
                                 // temperature T0
  // Returns adiabatic compression/expansion temperature corresponding
  // to given pressure
  var adiabaticTemp(const var& pres);
  
  
} // namespace tasystem



#endif // VARUTILS_H
//////////////////////////////////////////////////////////////////////

