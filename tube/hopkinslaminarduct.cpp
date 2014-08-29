#include "hopkinslaminarduct.h"
#include "tube.h"
#include "bessel.h"
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

    vd& xv=geom.xv;
    d& L=geom.L;
    vd Tmirror=Tl+(Tr-Tl)*math_common::skewsine(xv/L);
    // vd Tmirror=Tl+(Tr-Tl)*xv/L;
    se.setTMirror(Tmirror);
    d T;
    for(us i=0;i<vvertex.size();i++){
	TubeVertex& cvertex=*static_cast<TubeVertex*>(vvertex[i].get());
	T=Tmirror(i);
	cvertex.Ts.set(0,T);
	cvertex.T.set(0,T);
    }
    vd dTwdx=math_common::ddx_central(Tmirror,geom.xv);
    hopkinsheat.setdTwdx(geom,dTwdx);
  }
  
  
  
}
