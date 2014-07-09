#include "seg.h"
#include "bcvertex.h"

namespace tasystem{

  using segment::Seg;
  using segment::Vertex;
  using segment::BcVertex;
  using segment::connectpos;

  void connectbc(Seg&,const BcVertex&);

  Seg* copyseg(const Seg& orig);
  BcVertex* copybc(const BcVertex& orig);
} // namespace tasystem
