#include "hopkinslaminarduct.h"
#include "tube.h"

namespace tube{

  HopkinsLaminarDuct::HopkinsLaminarDuct(const HopkinsLaminarDuct& o):HopkinsLaminarDuct(o.geom)
  {
    LaminarDuct::operator=(o);
    hopkinsheat=o.hopkinsheat;    
  }
  HopkinsLaminarDuct& HopkinsLaminarDuct::operator=(const HopkinsLaminarDuct& o)
  {
    LaminarDuct::operator=(o);
    hopkinsheat=o.hopkinsheat;
    return *this;
  }
  HopkinsLaminarDuct HopkinsLaminarDuctTs(const Geom& geom,d Ts){
    HopkinsLaminarDuct d(geom);
    d.se.setTs(geom,Ts);
    return d;
  }
  SegBase* HopkinsLaminarDuct::copy() const{
    return new HopkinsLaminarDuct(*this);
  }

}
