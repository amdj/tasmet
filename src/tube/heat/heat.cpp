#include "heat.h"
#include "cell.h"
#include "jacrow.h"

namespace tube{

  vd HeatSource::Qsf(const Cell& v) const {return vd(v.gc->Ns(),fillwith::zeros);}
  tasystem::JacRow HeatSource::dQsf(const Cell& v) const {return tasystem::JacRow(-1,0);}

} // namespace tube

