#pragma once
#ifndef _HOPKINSHEAT_H_
#define _HOPKINSHEAT_H_
#include <vtypes.h>
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
    vc HeatTransferCoefH(const TubeVertex&) const;

    // This function computes for all nonzero frequencies the heat
    // transfer coeficient mathcalH.

  public:
    HopkinsHeatSource(const Tube& t);
    HopkinsHeatSource& operator=(const HopkinsHeatSource&);
    HopkinsHeatSource(const HopkinsHeatSource& o);

    virtual vd heat(const TubeVertex& v) const;
    // virtual dmat dUi(const TubeVertex& v) const;
    virtual dmat dTi(const TubeVertex& v) const;
  private:
    void setZeroFreq(const string&);
  };

} // namespace tube
#endif /* _HOPKINSHEAT_H_ */
