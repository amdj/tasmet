#include "tube.h"
#include "pressurebc.h"
#include "tasystem.h"
#include "globalconf.h"

namespace tube{
  using variable::var;

  using tasystem::TaSystem;
  using tasystem::Globalconf;

  PressureBc::PressureBc(const var& pres,const var& temp,us segnr,pos position):
    segnr(segnr),
    position(position),
    pLbc(pres),
    TLbc(temp)
  {
    TRACE(8,"PressureBc full constructor");
  }
  PressureBc::PressureBc(const var& pres,us segnr,pos position):
    segnr(segnr),
    position(position),
    pLbc(pres),
    TLbc(pres.gc())
  {
    TRACE(8,"PressureBc constructor for given pressure. Temperature computed");    
    const Globalconf& gc=pres.gc();
    d T0=gc.T0;
    d gamma=gc.gas.gamma(T0);
    vd p0(gc.Ns(),fillwith::ones); p0*=gc.p0;
    vd TLbct=T0*pow((p0+pLbc.tdata())/p0,(gamma-1.0)/gamma);		// Adiabatic compression/expansion
    TLbc.settdata(TLbct);
  }
  PressureBc::PressureBc(const PressureBc& other):
    PressureBc(other.pLbc,other.TLbc,other.segnr,other.position)
  {
    TRACE(8,"PressureBc copy constructor");
    TRACE(5,"pLbc:"<<pLbc());
  }
  void PressureBc::updateNf(){
    pLbc.updateNf();
    TLbc.updateNf();
  }

  void PressureBc::init(const TaSystem& sys)
  {
    TRACE(8,"PressureBc::initTubeVertex()");
    pLbc.setGc(sys.gc);
    TLbc.setGc(sys.gc);
  }


} // namespace tube



