// impedancebc.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////
// #define TRACERPLUS 20
#include "Python.h"
#include "impedancebc.h"
#include "weightfactors.h"

#include "duct.h"
#include "bccell.h"
#include "jacobian.h"
#include "pycallback.h"
#define fDFT (gc->fDFT)
#define iDFT (gc->iDFT)
#define DDTfd (gc->DDTfd)

#define Ns (gc->Ns())
#define Nf (gc->Nf())
#define eye (eye(Ns,Ns))

namespace duct {
  using namespace tasystem;
  
  // Return a var which contains the impedance using evaluation of
  // impfunc. It assumes a valid function, which should be tested
  // beforehand.
  var Zvar(PyObject* impfunc,const Globalconf* gc){
    TRACE(15,"Zvar()");
    vd omg=linspace(0,Nf*gc->getomg(),Nf+1);
    vc Zcomplex(Nf+1);
    VARTRACE(10,omg);
    for(us i=0;i<Nf+1;i++)
      Zcomplex(i)=python::PythonCallBackComplex(omg(i),impfunc);
    return var(*gc,Zcomplex);
  }
  void testImpedanceFunc(PyObject* impedancefunc){
    TRACE(15,"TestImpedanceFunc()");
    if(!impedancefunc)
      throw MyError("impedancefunc is null pointer!");
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
  }

  ImpedanceBc::ImpedanceBc(const string& segid,Pos position,PyObject* impedancefunc,d T0):
    DuctBc(segid,position),
    T0(T0)
  {
    TRACE(15,"ImpedanceBc::ImpedanceBc()");

    testImpedanceFunc(impedancefunc);
    // Increase reference to impedanceFunc, if not already error
    // thrown
    TRACE(15,"Incref impedance func");
    Py_INCREF(impedancefunc);
    impedanceFunc=impedancefunc;
  }
  ImpedanceBc::ImpedanceBc(const ImpedanceBc& other,const TaSystem& sys):
    DuctBc(other,sys),
    impedanceFunc(other.impedanceFunc),
    T0(other.T0)
  {
    TRACE(15,"ImpedanceBc::ImpedanceBc(copy and init)");
    assert(gc);
    testImpedanceFunc(impedanceFunc);
    // Increase reference to impedanceFunc
    TRACE(15,"Incref impedance func");
    Py_INCREF(impedanceFunc);    
  }
  ImpedanceBc::~ImpedanceBc(){
    TRACE(15,"ImpedanceBc::~ImpedanceBc()");
    // if(impedanceFunc)
      // Py_DECREF(impedanceFunc);
  }
  void ImpedanceBc::updateNf() {

  }
  void ImpedanceBc::setEqNrs(us firsteqnr) {
    TRACE(15,"void ImpedanceBc::setEqNrs()");
    VARTRACE(10,firsteqnr);
    DuctBc::setEqNrs(firsteqnr);
  } 

  vd ImpedanceBc::error() const {
    TRACE(15,"ImpedanceBc::error()");
    vd error(getNEqs());
    const BcCell& cell=t->bcCell(pos);

    // Impedance error
    TRACE(10,"Calling impedance function");
    dmat Z=Zvar(impedanceFunc,gc).freqMultiplyMat();
    TRACE(10,"Calling impedance function done.");
    // VARTRACE(15,Z);
    d Sf=cell.Sfbc();
    error.subvec(0,Ns-1)=cell.pbc()();
    error.subvec(0,Ns-1)+=-Z*cell.ubc()();

    // Isentropic state eq. error
    const vd& p=cell.pbc()();
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
    d Sf=pos==Pos::left?cell.Sfl:cell.Sfr;
    impjac+=JacCol(cell.pbc(),eye);
    impjac+=JacCol(cell.ubc(),-Z);
    jac+=impjac;


    // Adiabatic pressure-temperature Jacobian:
    JacRow IsentropicJac(firsteqnr+Ns,2);
    const var& Tbc=cell.Tbc();

    // Pressure part
    IsentropicJac+=JacCol(cell.pbc(),eye);
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
    const char* side;
    if(pos==Pos::left)
      side="left";
    else
      side="right";
    cout << "ImpedanceBc boundary condition set at "<<side <<" side of segment "<< segid<<".\n";
    if(detailnr>1){
      cout << "Prescribed impedance:" << Zvar(impedanceFunc,gc).getcRes() << "\n";
      cout << "Prescribed fluid temperature:" << T0 << "\n";      
      // cout << "Prescribed solid temperature:" << prescribeTs.getVals()() << "\n";      
    }
  } // show

} // namespace duct

//////////////////////////////////////////////////////////////////////
