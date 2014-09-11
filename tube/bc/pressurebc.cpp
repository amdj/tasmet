#include "tube.h"
#include "pressurebc.h"


// Can be removed later
#include "energyeq.h"
namespace tube{
  
  vd LeftPressureEq::error(const TubeVertex& v) const{
    TRACE(8,"LeftPressureEq::error()");
    return  v.pL()()-pLbc();
  }
  JacRow LeftPressureEq::jac(const TubeVertex& v) const{
    TRACE(8,"LeftPressureEq::jac()");
    JacRow jac(dofnr);
    jac+=dpL(v);
    return jac;
  }
  JacCol LeftPressureEq::dpL(const TubeVertex& v) const {
    TRACE(6,"LeftPressureEq::dpL()");
    return JacCol(v.pL(),arma::eye<dmat>(v.gc->Ns,v.gc->Ns));
  }
  LeftPressure::LeftPressure(const var& pres,const var& temp):
    TubeBcVertex(),
    pLbc(pres),
    TLbc(temp),
    leq(pLbc)
  {
    TRACE(8,"LeftPressure full constructor");
  }
  LeftPressure::LeftPressure(const var& pres):
    TubeBcVertex(),
    pLbc(pres),
    TLbc(*pres.gc),
    leq(pLbc)
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
  }
  void LeftPressure::initTubeVertex(us i,const Tube& thisseg)
  {
    TRACE(8,"LeftPressure::initTubeVertex()");
    TubeVertex::initTubeVertex(i,thisseg);
    pLbc.gc=thisseg.gc;
    TLbc.gc=thisseg.gc;
    eqs.push_back(std::unique_ptr<TubeEquation>(leq.copy()));
    // TRACE(100,"TLbc:"<<TLbc());
    // if(eqs.at(2)->getType()!=EqType::Ise){
    //   TRACE(100,"Changing energy equation...");
    //   // eq[2]=&peq;
    // }
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

    cWim1=0;
    cWi=w.wRl-w.wL0;
    cWip1=w.wRr-w.wL1;

    // Change momentum equation for open boundary, and prescribed pressure
    mWuim1=0;
    mWui=w.wRl/w.vSfR-w.wL0/w.vSfL;
    mWuip1=w.wRr/w.vSfR-w.wL1/w.vSfL;

    WARN("Not yet updated!");
    // mWpim1=0;     
    // mWpi=w.vSf*w.wRl;
    // mWpip1=w.vSf*w.wRr;
    // Change energy equation for open boundary and prescribed pressure
    eWgim1= 0;
    eWgi=   w.wRl-w.wL0;
    eWgip1= w.wRr-w.wL1;

    eWkinim1=0;
    eWkini=w.wRl/pow(w.vSfR,2)-w.wL0/pow(w.vSfL,2);    
    eWkinip1=w.wRr/pow(w.vSfR,2)-w.wL1/pow(w.vSfL,2);    

    TRACE(8,"Sf:"<<w.vSf);
    d x0=w.xvi;
    d x1=w.xvip1;
    TRACE(8,"x0:"<<x1);
    TRACE(8,"x1:"<<x1);
    d denom=x0*(1-x0/x1);
    TRACE(8,"denom:"<<denom);
    d x0_ov_x1sq=pow(x0/x1,2);
    eWc1=0;
    eWc2=lg.SfL/x0;
    eWc3=lg.vSf/w.dxp;
    eWc4=-w.vSfR/w.dxp;    

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
    d TLbc0=TLbc(0);
    // TRACE(100,"TLbc(0):"<<TLbc0);
    const Energy& e=static_cast<const Energy&>(*eqs[2]);
    d gamma=e.gamma(*this);
    d gamfac=gamma/(gamma-1.0);
    
    // TRACE(100,"TLbc:"<<TLbc());
    const LocalGeom rlg(lg.geom->localGeom(i+1));
    d xp2=rlg.xvi;
    d xp1=lg.xvi;

    vd T0=gc->T0*vd(gc->Ns,fillwith::ones);
    vd kappaL=e.kappaL(*this);

    d x0=lg.xvi;
    d x1=lg.xvi+w.dxp;
    d denom=x0*(1.0-x0/x1);

    d x0_ov_x1sq=pow(x0/x1,2);
    TRACE(12,"esource:"<<esource);
    esource+=-1.0*lg.SfL*fDFT*(kappaL%TLbct)/x0;
    TRACE(12,"esource:"<<esource);
    return esource;    
  }

} // namespace tube



