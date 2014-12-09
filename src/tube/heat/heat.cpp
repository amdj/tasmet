#include "tubevertex.h"

namespace tube{

  vd HeatSource::heat(const TubeVertex& v) const {return vd(v.gc->Ns(),fillwith::zeros);}
  dmat HeatSource::dUi(const TubeVertex& v) const {return zeros<dmat>(v.gc->Ns(),v.gc->Ns());}
  // dmat HeatSource::dpi(const TubeVertex& v) const {return v.zero;}  
  dmat HeatSource::dTi(const TubeVertex& v) const {return zeros<dmat>(v.gc->Ns(),v.gc->Ns());}



} // namespace tube

