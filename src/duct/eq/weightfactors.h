// weightfactors.h
//
// Author: J.A. de Jong 
//
// Description:
// Returns weight factors WlL,WlR,WrL,WrR in that order.
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef _WEIGHTFACTORS_H_
#define _WEIGHTFACTORS_H_
#include <tuple>
#include "vtypes.h"

namespace duct{
  SPOILNAMESPACE
  class Cell;


  class WeightFactors{
  public:
    d WrR=0,WrL=0,WlR=0,WlL=0;
    WeightFactors(const Cell& c);
    void show() const;
  };

  // weightfactors for interpolating to cell wall (left side of cell)
  std::tuple<d,d> weightFactorsWall(const Cell& c);

  // anonymous weightfactors for extrapolation from cell walls at
  // interior to cell wall at end of duct.
  std::tuple<d,d> BcWeightFactorsW(const Cell& c);

  // anonymous weightfactors for extrapolation from inner vertices to
  // cell wall
  std::tuple<d,d> BcWeightFactorsV(const Cell& c);
}
#endif /* _WEIGHTFACTORS_H_ */

//////////////////////////////////////////////////////////////////////
