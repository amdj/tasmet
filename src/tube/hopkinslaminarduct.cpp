#include "cell.h"
#include "hopkinslaminarduct.h"
#include "geom.h"
#include "solidenergy.h"


namespace tube{
  using tasystem::Globalconf;
  using tasystem::TaSystem;

  HopkinsLaminarDuct::HopkinsLaminarDuct(const Geom& geom,d Tl):
    HopkinsLaminarDuct(geom,Tl,Tl){
    TRACE(15,"HopkinsLaminarDuct::HopkinsLaminarDuct()");
  }

  HopkinsLaminarDuct::HopkinsLaminarDuct(const Geom& geom,d Tl,d Tr):
    LaminarDuct(geom),hopkinsheat(*this),
    Tl(Tl),
    Tr(Tr)
  {
    TRACE(15,"HopkinsLaminarDuct::HopkinsLaminarDuct(Geom,Tl,Tr)");
  }  
  HopkinsLaminarDuct::HopkinsLaminarDuct(const HopkinsLaminarDuct& o,
                                         const TaSystem& sys):
    LaminarDuct(o,sys),
    hopkinsheat(*this),
    Tl(o.Tl),
    Tr(o.Tr)
  {
    TRACE(15,"HopkinsLaminarDuct::HopkinsLaminarDuct(other)");

    // Set time-avg data to make solving bit easier
    assert(cells.size()>0);
    const vd& vx=geom().vx_vec();
    const d& L=geom().L();
    // Tmirror=Tl+(Tr-Tl)*math_common::skewsine(xv/L);
    vd Tmirror=Tl+(Tr-Tl)*vx/L;
    vd dTwdx=((Tr-Tl)/L)*ones(vx.size());
    // WARN("ToBECANGED!!")
    d T;
    for(us i=0;i<cells.size();i++){
      T=Tmirror(i);
      Cell& ccell=*cells[i];
      tasystem::var Tvar(*gc);
      Tvar.setadata(0,T);
      static_cast<SolidTPrescribed*>(ccell.Eq(Sol))->setTs(T);
      ccell.setResVar(Varnr::T,Tvar);
      ccell.setResVar(Varnr::Ts,Tvar);
    }
    hopkinsheat.setdTwdx(dTwdx);
    setInit(true);

  }
  segment::Seg* HopkinsLaminarDuct::copy(const TaSystem& sys) const{
    TRACE(20,"HopkinsLaminarDuct::copy()");
    return new HopkinsLaminarDuct(*this,sys);
  }  
  
}

