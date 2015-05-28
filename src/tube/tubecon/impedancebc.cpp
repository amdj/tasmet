// impedancebc.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////
// #define TRACERPLUS 20
#include "impedancebc.h"
#include "weightfactors.h"

#include "tube.h"
#include "bccell.h"
#include "jacobian.h"
#include "pycallback.h"
#include "Python.h"
#define fDFT (gc->fDFT)
#define iDFT (gc->iDFT)
#define DDTfd (gc->DDTfd)

#define Ns (gc->Ns())
#define Nf (gc->Nf())
#define eye (eye(Ns,Ns))

namespace tube {
  using namespace tasystem;
  
  var Zvar(PyObject* impfunc,const Globalconf* gc){

    vd omg=linspace(0,Nf*gc->getomg(),Nf+1);
    vc Zcomplex(Nf+1);
    VARTRACE(10,omg);
    for(us i=0;i<Nf+1;i++)
      Zcomplex(i)=python::PythonCallBackComplex(omg(i),impfunc);
    return var(*gc,Zcomplex);
  }


  ImpedanceBc::ImpedanceBc(us segnr,Pos position,PyObject* impedancefunc,d T0):
    TubeBc(segnr,position),
    T0(T0)
  {
    TRACE(15,"ImpedanceBc::ImpedanceBc()");
    // if(!impedanceFunc)
    // Simply testing the impedanceFunc
    try{ 
      c test=python::PythonCallBackComplex(100,impedancefunc);
      TRACE(15,"test:"<<test);
    }
    catch(MyError& m) {
      WARN("Something wrong with impedance callback function:");
      WARN(m.what());
      throw MyError("Illegal impedance callback function.");
    }
      
    // Increase reference to impedanceFunc, if not already error
    // thrown
    impedanceFunc=impedancefunc;
    Py_INCREF(impedanceFunc);
  }
  ImpedanceBc::ImpedanceBc(const ImpedanceBc& other,const TaSystem& sys):
    TubeBc(other,sys),
    impedanceFunc(other.impedanceFunc),
    T0(other.T0)
  {
    TRACE(15,"ImpedanceBc::ImpedanceBc(copy and init)");
    assert(gc);
    
    // Increase reference to impedanceFunc
    Py_INCREF(impedanceFunc);
    setInit(true);
  }
  ImpedanceBc::~ImpedanceBc(){
    if(impedanceFunc)
      Py_DECREF(impedanceFunc);
  }
  void ImpedanceBc::updateNf() {

  }
  void ImpedanceBc::setEqNrs(us firsteqnr) {
    TRACE(15,"void ImpedanceBc::setEqNrs()");
    VARTRACE(10,firsteqnr);
    TubeBc::setEqNrs(firsteqnr);
  } 

  vd ImpedanceBc::error() const {
    TRACE(15,"ImpedanceBc::error()");
    vd error(getNEqs());
    const BcCell& cell=t->bcCell(pos);

    // Impedance error
    dmat Z=Zvar(impedanceFunc,gc).freqMultiplyMat();
    // VARTRACE(15,Z);
    d Sf=pos==Pos::left?cell.SfL:cell.SfR;
    error.subvec(0,Ns-1)=cell.extrapolateQuant(Varnr::p);
    error.subvec(0,Ns-1)+=-Z*fDFT*(cell.mbc().tdata()/(Sf*iDFT*cell.extrapolateQuant(Varnr::rho)));

    // Isentropic state eq. error
    const vd p=cell.extrapolateQuant(Varnr::p);
    const var& Tbc=cell.Tbc();
    const d p0=gc->p0();
    const d gamma=gc->gas().gamma(T0);
    vd p0t(Ns,fillwith::ones); p0t*=p0;
    // Adiabatic compression/expansion
    error.subvec(Ns,2*Ns-1)=p-fDFT*(p0*(pow(Tbc.tdata()/T0,gamma/(gamma-1))-1));

    // Extrapolated enthalpy flow
    error.subvec(2*Ns,3*Ns-1)=-cell.mHbc()();
    error.subvec(2*Ns,3*Ns-1)+=cell.extrapolateQuant(Varnr::mH);
    return error;
  }
  void ImpedanceBc::jac(Jacobian& jac) const{
    TRACE(8,"ImpedanceBc::jac()");
    VARTRACE(10,firsteqnr);
    const BcCell& cell=t->bcCell(pos);

    // Impedance Jacobian
    // VARTRACE(25,cell.extrapolateQuant(Varnr::rho));
    dmat Z=Zvar(impedanceFunc,gc).freqMultiplyMat();
    JacRow impjac(firsteqnr,5); // 5=(2*p,2*rho,1*mbc)
    d Sf=pos==Pos::left?cell.SfL:cell.SfR;
    impjac+=cell.dExtrapolateQuant(Varnr::p);
    const vd rhobct=iDFT*cell.extrapolateQuant(Varnr::rho);
    impjac+=JacCol(cell.mbc(),-Z*fDFT*diagmat(1/(Sf*rhobct))*iDFT);
    d w1,w2; std::tie(w1,w2)=BcWeightFactorsV(cell);
    if(pos==Pos::left) {
      impjac+=JacCol(cell.rho(),w1*Z*fDFT*diagmat(cell.mbc().tdata()/(Sf*pow(rhobct,2)))*iDFT);      
      impjac+=JacCol(cell.rhoR(),w2*Z*fDFT*diagmat(cell.mbc().tdata()/(Sf*pow(rhobct,2)))*iDFT);      
    } else {
      impjac+=JacCol(cell.rho(),w1*Z*fDFT*diagmat(cell.mbc().tdata()/(Sf*pow(rhobct,2)))*iDFT);      
      impjac+=JacCol(cell.rhoL(),w2*Z*fDFT*diagmat(cell.mbc().tdata()/(Sf*pow(rhobct,2)))*iDFT);      
    }

    jac+=impjac;


    // Adiabatic pressure-temperature Jacobian:
    JacRow IsentropicJac(firsteqnr+Ns,2);
    const var& Tbc=cell.Tbc();

    // Pressure part
    IsentropicJac+=cell.dExtrapolateQuant(Varnr::p);
    const d p0=gc->p0();
    const d gamma=gc->gas().gamma(T0);

    // Temperature part
    d fac1=p0*gamma/(T0*(gamma-1));
    IsentropicJac+=JacCol(cell.Tbc(),-fDFT*diagmat(fac1*pow(Tbc.tdata()/T0,1/(gamma-1)))*iDFT);
    jac+=IsentropicJac;

    // Prescribed enthalpy flow.
    JacRow enthalpy_extrapolated_jac(firsteqnr+2*Ns,3);
    enthalpy_extrapolated_jac+=cell.dExtrapolateQuant(Varnr::mH);
    enthalpy_extrapolated_jac+=JacCol(cell.mHbc(),-eye);
    jac+=enthalpy_extrapolated_jac;

  }
  void ImpedanceBc::show(us detailnr) const {
    TRACE(5,"ImpedanceBc::show()");
    checkInit();
    const char* side;
    if(pos==Pos::left)
      side="left";
    else
      side="right";
    cout << "ImpedanceBc boundary condition set at "<<side <<" side of segment "<< segnr<<".\n";
    if(detailnr>1){
      cout << "Prescribed impedance:" << Zvar(impedanceFunc,gc).getcRes() << "\n";
      cout << "Prescribed fluid temperature:" << T0 << "\n";      
      // cout << "Prescribed solid temperature:" << prescribeTs.getVals()() << "\n";      
    }
  } // show

} // namespace tube

//////////////////////////////////////////////////////////////////////
