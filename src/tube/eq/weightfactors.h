// weightfactors.h
//
// Author: J.A. de Jong 
//
// Description:
// Returns weight factors WLl,WLr,WRl,WRr in that order.
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef _WEIGHTFACTORS_H_
#define _WEIGHTFACTORS_H_
#include <tuple>
#include "vtypes.h"

namespace tube{
  SPOILNAMESPACE
  class Cell;
  std::tuple<d,d,d,d> WeightFactors(const Cell& c);
}
#endif /* _WEIGHTFACTORS_H_ */

//////////////////////////////////////////////////////////////////////
