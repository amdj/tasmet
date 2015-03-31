#include "cell.h"

namespace tube{

  vd HeatSource::heat(const Cell& v) const {return vd(v.gc->Ns(),fillwith::zeros);}
  dmat HeatSource::dUi(const Cell& v) const {return zeros<dmat>(v.gc->Ns(),v.gc->Ns());}
  // dmat HeatSource::dpi(const Cell& v) const {return v.zero;}  
  dmat HeatSource::dTi(const Cell& v) const {return zeros<dmat>(v.gc->Ns(),v.gc->Ns());}



} // namespace tube

