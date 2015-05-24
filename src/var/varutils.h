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
                                 // temperature gc->T0


  // Returns adiabatic compression/expansion temperature corresponding
  // to given pressure. T0 is Grabbed from the Global configuration
  // (Globalconf) if an invalid value is given. The default argument
  // is an invalid value.
  var adiabaticTemp(const var& pres,d T0=-1);
  
  
} // namespace tasystem



#endif // VARUTILS_H
//////////////////////////////////////////////////////////////////////

