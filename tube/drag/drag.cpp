#include "tube.h"
#include "drag.h"
#include "tubevertex.h"
#include "segbase.h"
namespace tube{

  vd DragResistance::drag(const TubeVertex& v) const {return vd(v.gc->Ns(),fillwith::zeros);}
  dmat DragResistance::dUi(const TubeVertex& v) const {return zeros<dmat>(v.gc->Ns(),v.gc->Ns());}
  // dmat DragResistance::drhoi(const TubeVertex& v) const {return v.zero;}
  // dmat DragResistance::dpi(const TubeVertex& v) const {return v.zero;}


} // namespace tube





















