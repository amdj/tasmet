#include "hopkinslaminarduct.h"
#include "tube.h"

namespace tube{


  HopkinsLaminarDuct ColdHopkinsLaminarDuct(const Geom& geom,const Globalconf& gc){
    HopkinsLaminarDuct d(geom);
     d.se.setTs(gc.T0);
    return d;
  }
  

}
