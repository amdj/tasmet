#pragma once
#ifndef _HOPKINSHEAT_H_
#define _HOPKINSHEAT_H_
#include "vtypes.h"
#include "heat.h"

#include "rottfuncs.h"




namespace tube{
  SPOILNAMESPACE
  class Cell;
  class Tube;
  class Geom;
  
  class HopkinsHeatSource:public HeatSource{
    string cshape;
    d (*zeroheatH_funptr)(d,d)=nullptr;
    d zeroheatQ=0;
    rottfuncs::RottFuncs rf;
    // This function computes for all nonzero frequencies the heat
    // transfer coeficient mathcalH.
    vd dTwdx;
  public:
    HopkinsHeatSource(const Tube& t);
    HopkinsHeatSource& operator=(const HopkinsHeatSource&)=delete;
    HopkinsHeatSource(const HopkinsHeatSource& o)=delete;
    void setdTwdx(const vd& dTwdx);
    vd heat(const Cell& v) const;
    dmat dmi(const Cell& v) const;
    dmat dTi(const Cell& v) const;
  // private:
    void setZeroFreq(const string&);
    vc HeatTransferCoefH(const Cell&) const;
    vc HeatTransferCoefQ(const Cell&) const;    
  };

} // namespace tube
#endif /* _HOPKINSHEAT_H_ */
