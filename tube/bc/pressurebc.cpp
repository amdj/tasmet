#include "tube.h"
#include "pressurebc.h"

// Can be removed later
#include "energyeq.h"
namespace tube{
  
  JacCol pl(const TubeVertex& v){
    JacCol dpL(v.pL());
    dpL.setToAdd(false);        // We do not want dpL() anymore
    return dpL;
  }
  JacCol LeftPressureMomentumEq::dpL(const TubeVertex& v) const {
    TRACE(6,"LeftPressureMomentumEq::dpL()");
    return pl(v);
  }
  JacCol LeftPressureEnergyEq::dpL(const TubeVertex& v) const {
    TRACE(6,"LeftPressureEnergyEq::dpL()");
    return pl(v);
  }
  JacCol LeftPressureIsentropicEq::dpL(const TubeVertex& v) const {
    TRACE(6,"LeftPressureIsentropicEq::dpL()");
    return pl(v);
  }
  JacCol LeftPressureStateEq::dpL(const TubeVertex& v) const {
    TRACE(6,"LeftPressureStateEq::dpL()");
    return pl(v);
  }
  LeftPressure::LeftPressure(const var& pres,const var& temp):
    pLbc(pres),
    TLbc(temp)
  {
    TRACE(8,"LeftPressure full constructor");
  }
  LeftPressure::LeftPressure(const var& pres):
    pLbc(pres),
    TLbc(*pres.gc)
  {
    TRACE(8,"LeftPressure constructor for given pressure. Temperature computed");    
    const Globalconf* gc=pres.gc;
    d T0=gc->T0;
    d gamma=gc->gas.gamma(T0);
    vd p0(gc->Ns,fillwith::ones); p0*=gc->p0;
    vd TLbct=T0*pow((p0+pLbc.tdata())/p0,(gamma-1.0)/gamma);		// Adiabatic compression/expansion
    TLbc.settdata(TLbct);
    // TRACE(100,"TLbc(0):"<<TLbc);
  }
  LeftPressure::LeftPressure(const LeftPressure& other):
    LeftPressure(other.pLbc,other.TLbc)
  {
    TRACE(8,"LeftPressure copy constructor");
    TRACE(15,"pLbc:"<<pLbc());
  }
  void LeftPressure::initTubeVertex(us i,const Tube& thisseg)
  {
    TRACE(8,"LeftPressure::initTubeVertex()");
    TubeVertex::initTubeVertex(i,thisseg);
    pLbc.gc=thisseg.gc;
    TLbc.gc=thisseg.gc;
    // eqs.at(3).reset(leq.copy()); // Replace equation of state for the
    // boundary condition on pressure
    lmomeq.init(thisseg);
    leneq.init(thisseg);
    lseq.init(thisseg);

    // Set our own momentumeq,energy eq and state eq
    eqs.at(1)=&lmomeq;
    if(thisseg.getName().compare("IsentropicTube")==0){
      TRACE(15,"Tube is isentropic, adapting isentropic state for leftpressure bc");
      eqs.at(2)=&liseq;
    }
    else{
      eqs.at(2)=&leneq;    // Remove state equation from list
      TRACE(15,"Tube is not isentropic, full energy eq for leftpressure bc");
    }
    eqs.erase(eqs.begin()+3);

    // Refill vars
    vars.clear();
    vars.push_back(&rho);
    vars.push_back(&U);
    // vars.push_back(&p); // No p here
    vars.push_back(&T);
    vars.push_back(&Ts);    
    assert(vars.size()==4);
    
    LeftPressure::updateW(thisseg);
  }
  void LeftPressure::show() const {
    cout << "LeftPressure boundary condition set at left side of tube.\n";
    cout << "Number of equations: " << eqs.size() << "\n";
    TubeVertex::show();
  }
  void LeftPressure::updateW(const SegBase& thisseg)
  {
    // Change continuity equation for an open boundary
    TRACE(8,"LeftPressure::updateW()");

    c.Wim1=0;
    c.Wi=wRl-wL0;
    c.Wip1=wRr-wL1;

    const LocalGeom& rlg=right->lg;
    
    // Change momentum equation for open boundary, and prescribed pressure

    lmomeq.Wddt=m.Wddt;
    lmomeq.Wuim1=0;
    lmomeq.Wui=wRl/lg.vSf-wL0/lg.SfL;
    lmomeq.Wuip1=wRr/rlg.vSf-wL1/lg.SfL;
    lmomeq.WpL=m.WpL;
    lmomeq.WpR=m.WpR;

    leneq.Wddt=e.Wddt;
    leneq.Wddtkin=e.Wddtkin;    
    // Change energy equation for open boundary and prescribed pressure
    leneq.Wgim1= 0;
    // Left boundary, velocity 
    leneq.Wgim=-wL0;
    leneq.WgUip1pL=-wL1;
    
    leneq.Wgip=wRl;
    leneq.Wgip1= wRr;    

    leneq.Wkinim1=0;
    leneq.Wkini=wRl/pow(lg.vSf,2)-wL0/pow(lg.SfL,2);    
    leneq.Wkinip1=wRr/pow(rlg.vSf,2)-wL1/pow(lg.SfL,2);    

    d x0=xvi;

    leneq.Wc1=0;
    leneq.Wc2=lg.SfL/x0;
    leneq.Wc3=lg.SfR/dxp;
    leneq.Wc4=-lg.SfR/dxp;    

  }
  vd LeftPressure::msource() const{
    TRACE(5,"LeftPressure::msource()");
    vd msource(gc->Ns,fillwith::zeros);
    // msource=-1.0*lg.SfL*pLbc();
    return msource;
  }
  vd LeftPressure::esource() const {
    TRACE(5,"LeftPressure::esource()");
    vd esource=TubeVertex::esource();
    const dmat& fDFT=gc->fDFT;
    // TRACE(100,"Stupid hack to test energy source");
    vd TLbct=TLbc.tdata();

    const Energy& e=static_cast<const Energy&>(*eqs[2]);
    
    vd kappaL=leneq.kappaL(*this);

    d x0=lg.xvi;

    esource+=-1.0*lg.SfL*fDFT*(kappaL%TLbct)/x0;
    TRACE(12,"esource:"<<esource);
    return esource;    
  }

} // namespace tube



