// hopkinsheat.h
//
// Author: J.A. de Jong 
//
// Description:
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef HOPKINSHEAT_H
#define HOPKINSHEAT_H

#include "vtypes.h"
#include "heat.h"
#include "rottfuncs.h"

namespace tasystem {
  class JacRow;
} // namespace tasystem
namespace tube{
  SPOILNAMESPACE
  class Cell;
  class LaminarDuct;
  class Geom;
  
  class HopkinsHeatSource:public HeatSource{
    string cshape;
    //  Pointer to function that computes heat transfer from wall to
    //  fluid at frequency zero
    d (*zeroheatH_funptr)(d,d)=nullptr;
    d zeroheatQ=0;
    const LaminarDuct* t;
    rottfuncs::RottFuncs rf;
  public:
    HopkinsHeatSource(const LaminarDuct& t);
    HopkinsHeatSource& operator=(const HopkinsHeatSource&)=delete;
    HopkinsHeatSource(const HopkinsHeatSource& o)=delete;
    vd Qsf(const Cell& v) const;
    tasystem::JacRow dQsf(const Cell&) const;
  private:
    void setZeroFreq(const string&);
    vc HeatTransferCoefH(const Cell&) const;
    vc HeatTransferCoefQ(const Cell&) const;
    d dTwdx(const Cell&) const;
  };

} // namespace tube

#endif // HOPKINSHEAT_H
//////////////////////////////////////////////////////////////////////
