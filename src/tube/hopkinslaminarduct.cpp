#include "tube.h"
#include "hopkinslaminarduct.h"
#include "bessel.h"
#include "geom.h"
#include "cell.h"
#include "constants.h"

namespace tube{
  using tasystem::Globalconf;
  using tasystem::TaSystem;

  HopkinsLaminarDuct::HopkinsLaminarDuct(const Geom& geom,d Tl):
    HopkinsLaminarDuct(geom,Tl,Tl){
    TRACE(15,"HopkinsLaminarDuct::HopkinsLaminarDuct()");
  }

  HopkinsLaminarDuct::HopkinsLaminarDuct(const Geom& geom,d Tl,d Tr):
    LaminarDuct(geom),hopkinsheat(*this){
    TRACE(15,"HopkinsLaminarDuct::HopkinsLaminarDuct(Geom,Tl,Tr)");
    this->Tl=Tl;
    this->Tr=Tr;
  }  
  HopkinsLaminarDuct::HopkinsLaminarDuct(const HopkinsLaminarDuct& o):LaminarDuct(o),
                                                                      hopkinsheat(o.hopkinsheat),
                                                                      Tl(o.Tl),
                                                                      Tr(o.Tr)
  {
    TRACE(15,"HopkinsLaminarDuct::HopkinsLaminarDuct(other)");

  }
  
  void HopkinsLaminarDuct::init(const TaSystem& sys){
    TRACE(15,"HopkinsLaminarDuct::init(gc)");
    LaminarDuct::init(sys);

    // Set time-avg data to make solving bit easier
    assert(cells.size()>0);
    vd vx(geom().nCells());
    for(us i=0;i<vx.size();i++)
      vx(i)=geom().vx(i);
    const d& L=geom().L();
    // Tmirror=Tl+(Tr-Tl)*math_common::skewsine(xv/L);
    Tmirror=Tl+(Tr-Tl)*vx/L;
    dTwdx=math_common::ddx_central(Tmirror,vx);
    // WARN("ToBECANGED!!")
    d T;
    for(us i=0;i<cells.size();i++){
      T=Tmirror(i);
      Cell& ccell=*cells[i];
      variable::var Tvar(*gc);
      Tvar.set(0,T);
      ccell.setResVar(varnr::T,Tvar);
      ccell.setResVar(varnr::Ts,Tvar);
    }
    hopkinsheat.setdTwdx(geom(),dTwdx);
    setInit(true);
  }
  
}

