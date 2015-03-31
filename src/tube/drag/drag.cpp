#include "drag.h"
#include "cell.h"

namespace tube{

  vd DragResistance::drag(const Cell& v) const {return vd(v.gc->Ns(),fillwith::zeros);}
  dmat DragResistance::dUi(const Cell& v) const {return zeros<dmat>(v.gc->Ns(),v.gc->Ns());}
  // dmat DragResistance::drhoi(const Cell& v) const {return v.zero;}
  // dmat DragResistance::dpi(const Cell& v) const {return v.zero;}


} // namespace tube





















