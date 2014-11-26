#include "tube.h"
#include "hopkinslaminarduct.h"
#include "bessel.h"
#include "tubevertex.h"

namespace tube{
  HopkinsLaminarDuct::HopkinsLaminarDuct(const Geom& geom,d Tl):HopkinsLaminarDuct(geom,Tl,Tl){}

  HopkinsLaminarDuct::HopkinsLaminarDuct(const Geom& geom,d Tl,d Tr):LaminarDuct(geom),hopkinsheat(*this){
    this->Tl=Tl;
    this->Tr=Tr;
  }  
  HopkinsLaminarDuct::HopkinsLaminarDuct(const HopkinsLaminarDuct& o):LaminarDuct(o),
										hopkinsheat(o.hopkinsheat),
										Tl(o.Tl),
										Tr(o.Tr)
  {  }
  HopkinsLaminarDuct& HopkinsLaminarDuct::operator=(const HopkinsLaminarDuct& o)
  {
    LaminarDuct::operator=(o);
    Tl=o.Tl;
    Tr=o.Tr;
    hopkinsheat=o.hopkinsheat;
    
    return *this;
  }
  
  void HopkinsLaminarDuct::init(const Globalconf& gc){
    TRACE(15,"HopkinsLaminarDuct::init(gc)");
    LaminarDuct::init(gc);
    // Set time-avg data to make solving bit easier
    assert(vvertex.size()>0);
    vd vx(geom().nCells());
    for(us i=0;i<vx.size();i++)
      vx(i)=geom().vx(i);
    const d& L=geom().L();
    // Tmirror=Tl+(Tr-Tl)*math_common::skewsine(xv/L);
    Tmirror=Tl+(Tr-Tl)*vx/L;
    dTwdx=math_common::ddx_central(Tmirror,vx);
    // WARN("ToBECANGED!!")
    d T;
    for(us i=0;i<vvertex.size();i++){
      T=Tmirror(i);
      TubeVertex& cvertex=*static_cast<TubeVertex*>(vvertex[i]);
      cvertex.Ts.set(0,T);
      cvertex.T.set(0,T);
    }
    hopkinsheat.setdTwdx(geom(),dTwdx);
  }
  
}

