// #define TRACERPLUS 20

#include "tube.h"
#include "bccell.h"
#include "weightfactors.h"
#include "jacrow.h"
#include "var.h"
#include "momentum.h"

#include "drag.h"

#ifdef NODRAG
#warning Drag is turned off!
#endif

#define iDFT (v.gc->iDFT)
#define fDFT (v.gc->fDFT)
#define DDTfd (v.gc->DDTfd)
#define Ns (v.gc->Ns())
#define eye (arma::eye(Ns,Ns))

namespace tube{
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;
  using tasystem::var;


  void Momentum::show() const{
    cout << "----------------- Momentum equation\n";
    cout << "Wddt   : " << Wddt << "\n";
    cout << "Wpi    : " << Wpi << "\n";
    cout << "Wpm1    : " << Wpim1 << "\n";    
  }
  Momentum::~Momentum(){TRACE(0,"~Momentum()");}
  void Momentum::init()  {
    TRACE(5,"Momentum::init(tube)");
    t=&v.getTube();
    if(v.left()){
      Wddt=v.vx-v.left()->vx;;
      Wpi=v.Sfl;
      Wpim1=-v.Sfl;
    }
  }
  vd Momentum::error() const {		// Error in momentum equation
    TRACE(15,"Momentum::Error()");

    // Solve momentum equation only for interior walls
    assert(v.i>0);
    
    vd error=zeros();
    const vd& rhoti=v.rho().tdata();

    error+=Wddt*DDTfd*v.ml()();

    // Pressure right
    error+=Wpi*v.p()();
    // Pressure left
    error+=Wpim1*v.left()->p()();

    error+=v.mu()();
    error-=v.left()->mu()();
    // (m)^2 in time domain at i

    // Drag term
    #ifndef NODRAG
    assert(t);
    error+=Wddt*t->dragResistance().drag(v);
    #endif
    // (Boundary) source term
    error+=v.msource();
    return error;
  }
  JacRow Momentum::jac() const {
    TRACE(15,"Momentum::jac()");
    JacRow jac(dofnr,9);

    // Solve momentum equation only for interior walls
    assert(v.i>0);

    // Time-derivative of mass flow
    jac+=JacCol(v.ml(),Wddt*DDTfd);

    jac+=JacCol(v.p(),Wpi*eye);
    jac+=JacCol(v.left()->p(),Wpim1*eye);
    jac+=JacCol(v.mu(),eye);
    jac+=JacCol(v.left()->mu(),-eye);

 
    #ifndef NODRAG
    jac+=JacCol(v.ml(),Wddt*t->dragResistance().dm(v));
    #endif
   
    return jac;
  }
  void Momentum::domg(vd & domg_) const {
    TRACE(18,"Momentum::domg()");
    // Possibly later adding drag->domg();
    vd domg_full=(Wddt/v.gc->getomg())*DDTfd*v.ml()();
    domg_.subvec(dofnr,dofnr+Ns-1)=domg_full;
    TRACE(0,"Momentum::domg() done");
  }


  namespace LEFT {
  
   static vd2 weightfactors(const Cell& v){
      d vxi=v.vx;
      d vxip1=v.right()->vx;

      // Compute weight factors
      d wL0=vxip1/(vxip1-vxi);
      d wL1=-vxi/(vxip1-vxi);
      // VARTRACE(25,wL0);
      // VARTRACE(25,wL1);

      return vd2({wL0,wL1});
    }
  
    static vd extrapolatePressure(const Cell& v){
      TRACE(15,"LEFT::extrapolatePressure()");
      vd2 w=weightfactors(v); d wL0=w(0),wL1=w(1);
      assert(v.right());
      return wL0*v.p()()+wL1*v.right()->p()();
    }
    static vd extrapolateMomentumFlow(const Cell& v){
      vd2 w=weightfactors(v); d wL0=w(0),wL1=w(1);
      return wL0*v.mu()()+wL1*v.right()->mu()();
    }
    static JacRow dExtrapolateMomentumFlow(const Cell& v){
      vd2 w=weightfactors(v); d wL0=w(0),wL1=w(1);
      JacRow jacrow(2);
      jacrow+=JacCol(v.mu(),wL0*eye);
      jacrow+=JacCol(v.right()->mu(),wL1*eye);
      return jacrow;
    }
    static JacRow dExtrapolatePressure(const Cell& v){
      vd2 w=weightfactors(v); d wL0=w(0),wL1=w(1);
      JacRow jacrow(2);
      jacrow+=JacCol(v.p(),wL0*eye);
      jacrow+=JacCol(v.right()->p(),wL1*eye);
      return jacrow;
    }
  } // namespace LEFT

  namespace RIGHT {
   static std::tuple<d,d> weightfactors(const Cell& v){
      d xr=v.xr;
      d vxm1=v.vx;
      d vxm2=v.left()->vx;

      // Compute weight factors
      d wRNm1=(vxm2-xr)/(vxm2-vxm1);
      d wRNm2=(xr-vxm1)/(vxm2-vxm1);

      // VARTRACE(25,wRNm1);
      // VARTRACE(25,wRNm2);
      return std::make_tuple(wRNm1,wRNm2);
    }

   // Extrapolate momentum flow left side  
   static vd extrapolateMomentumFlow(const Cell& v){
      d wRNm1,wRNm2; std::tie(wRNm1,wRNm2)=weightfactors(v);
      assert(v.left());
      return wRNm1*v.mu()()+wRNm2*v.left()->mu()();
    }
   static vd extrapolatePressure(const Cell& v){
      TRACE(15,"RIGHT::extrapolatePressure()");
      d wRNm1,wRNm2; std::tie(wRNm1,wRNm2)=weightfactors(v);
      assert(v.left());
      return wRNm1*v.p()()+wRNm2*v.left()->p()();
    }
   static JacRow dExtrapolateMomentumFlow(const Cell& v){
      d wRNm1,wRNm2; std::tie(wRNm1,wRNm2)=weightfactors(v);
      assert(v.left());
      JacRow jacrow(2);
      jacrow+=JacCol(v.mu(),wRNm1*eye);
      jacrow+=JacCol(v.left()->mu(),wRNm2*eye);
      return jacrow;
    }

   static JacRow dExtrapolatePressure(const Cell& v){
      d wRNm1,wRNm2; std::tie(wRNm1,wRNm2)=weightfactors(v);
      assert(v.left());
      JacRow jacrow(2);
      jacrow+=JacCol(v.p(),wRNm1*eye);
      jacrow+=JacCol(v.left()->p(),wRNm2*eye);
      return jacrow;
    }

  } // namespace RIGHT

  ExtrapolatePressure::ExtrapolatePressure(const BcCell& v):Equation(v),v(v){}
  BcVelocity::BcVelocity(const BcCell& v):Equation(v),v(v){}

  vd ExtrapolatePressure::extrapolatePressure(const Cell& v) {
    TRACE(15,"ExtrapolatePressure::extrapolatePressure()");

    if(!v.left() && v.right())
      return LEFT::extrapolatePressure(v);
    else
      return RIGHT::extrapolatePressure(v);
  }
  vd ExtrapolateMomentumFlow::extrapolateMomentumFlow(const Cell& v) {
    TRACE(15,"ExtrapolateMomentumFlow::extrapolateMomentumFlow()");

    if(!v.left() && v.right())
      return LEFT::extrapolateMomentumFlow(v);
    else
      return RIGHT::extrapolateMomentumFlow(v);
  }
  vd ExtrapolatePressure::error() const {
    TRACE(15,"ExtrapolatePressure::error()");
    assert((!v.left() && v.right()) || (!v.right() && v.left()));
    VARTRACE(15,v.pbc()());
    return extrapolatePressure(v)-v.pbc()();
  }
  vd BcVelocity::error() const {
    TRACE(15,"BcVelocity::error()");
    assert((!v.left() && v.right()) || (!v.right() && v.left()));
   const vd& mbct=v.mbc().tdata();
    const vd& rhobct=v.rhobc().tdata();
    return fDFT*(mbct/(rhobct*v.Sfbc()))-v.ubc()();
  }
  JacRow ExtrapolatePressure::dExtrapolatePressure(const Cell& v) {
    TRACE(15,"ExtrapolatePressure::dExtrapolatePressure()");
    assert((!v.left() && v.right()) || (!v.right() && v.left()));
    if(!v.left() && v.right())
      return LEFT::dExtrapolatePressure(v);
    else
      return RIGHT::dExtrapolatePressure(v);
  }
  JacRow ExtrapolateMomentumFlow::dExtrapolateMomentumFlow(const Cell& v) {
    TRACE(15,"ExtrapolateMomentumFlow::dExtrapolateMomentumFlow()");
    assert((!v.left() && v.right()) || (!v.right() && v.left()));
    if(!v.left() && v.right())
      return LEFT::dExtrapolateMomentumFlow(v);
    else
      return RIGHT::dExtrapolateMomentumFlow(v);
  }

  JacRow ExtrapolatePressure::jac() const{
    TRACE(15,"ExtrapolatePressure::jac()");
    assert((!v.left() && v.right()) || (!v.right() && v.left()));
    JacRow dExtrapolatePressure(dofnr,3);
    dExtrapolatePressure+=this->dExtrapolatePressure(v);
    dExtrapolatePressure+=JacCol(v.pbc(),-eye);
    return dExtrapolatePressure;
  }
  JacRow BcVelocity::jac() const{
    TRACE(15,"BcVelocity::jac()");
    assert((!v.left() && v.right()) || (!v.right() && v.left()));
    JacRow dBcVelocity(dofnr,3);
    const vd& mbct=v.mbc().tdata();
    const vd& rhobct=v.rhobc().tdata();
    d Sfbc=v.Sfbc();
    dBcVelocity+=JacCol(v.mbc(),fDFT*diagmat(1/(rhobct*Sfbc))*iDFT);
    dBcVelocity+=JacCol(v.rhobc(),-fDFT*diagmat(mbct/(pow(rhobct,2))*Sfbc)*iDFT);
    dBcVelocity+=JacCol(v.ubc(),-eye);
    return dBcVelocity;
  }
  void ExtrapolatePressure::show() const{
    cout << "----------------- Extrapolate Pressure\n";
  }
 void ExtrapolatePressure::init() {
    TRACE(15,"ExtrapolatePressure::init()");
  }
  void BcVelocity::show() const{
    cout << "----------------- BcVelocity\n";
  }
 void BcVelocity::init() {
    TRACE(15,"BcVelocity::init()");
  }
  
} // namespace tube
