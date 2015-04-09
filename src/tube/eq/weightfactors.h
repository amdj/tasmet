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
  std::tuple<d,d> BcWeightFactors(const Cell& c);
}
#endif /* _WEIGHTFACTORS_H_ */

//////////////////////////////////////////////////////////////////////
