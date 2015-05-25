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


  class WeightFactors{
  public:
    d WRr=0,WRl=0,WLr=0,WLl=0;
    WeightFactors(const Cell& c);
  };
  // anonymous weightfactors for extrapolation from cell walls at
  // interior to cell wall at end of tube.
  std::tuple<d,d> BcWeightFactorsW(const Cell& c);

  // anonymous weightfactors for extrapolation from inner vertices to
  // cell wall
  std::tuple<d,d> BcWeightFactorsV(const Cell& c);
}
#endif /* _WEIGHTFACTORS_H_ */

//////////////////////////////////////////////////////////////////////
