#pragma once
#ifndef _HOPKINSHEAT_H_
#define _HOPKINSHEAT_H_
#include "vtypes.h"
#include "heat.h"
#include "tube.h"
#include "rottfuncs.h"

namespace tube{
  SPOILNAMESPACE
  class TubeVertex;
  // Stub for a DragResistance class



  
  class HopkinsHeatSource:public HeatSource{
    string cshape;
    d (*zeroheatH_funptr)(d,d);
    rottfuncs::RottFuncs rf;
    // This function computes for all nonzero frequencies the heat
    // transfer coeficient mathcalH.
    const vd* dTwdx;
  public:
    HopkinsHeatSource(const Tube& t);
    HopkinsHeatSource& operator=(const HopkinsHeatSource&);
    HopkinsHeatSource(const HopkinsHeatSource& o);
    void setdTwdx(const Geom& g,const vd& dTwdx);
    virtual vd heat(const TubeVertex& v) const;
    virtual dmat dUi(const TubeVertex& v) const;
    virtual dmat dTi(const TubeVertex& v) const;
  private:
    void setZeroFreq(const string&);
    vc HeatTransferCoefH(const TubeVertex&) const;
    vc HeatTransferCoefQ(const TubeVertex&) const;    
  };

} // namespace tube
#endif /* _HOPKINSHEAT_H_ */
