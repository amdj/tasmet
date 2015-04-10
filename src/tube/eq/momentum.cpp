// #define TRACERPLUS 20

#include "tube.h"
#include "cell.h"
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

namespace tube{
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;
  using variable::var;


  void Momentum::show() const{
    cout << "----------------- Momentum equation\n";
    cout << "Wddt   : " << Wddt << "\n";
    cout << "Wpi    : " << Wpi << "\n";
    cout << "Wpm1    : " << Wpim1 << "\n";    
  }
  Momentum::~Momentum(){TRACE(0,"~Momentum()");}
  void Momentum::init()
  {
    TRACE(5,"Momentum::init(tube)");
    t=&v.getTube();
    if(v.left()){
      Wddt=v.vx-v.left()->vx;;
      Wpi=v.SfL;
      Wpim1=-v.SfL;
    }
  }
  vd Momentum::error() const {		// Error in momentum equation
    TRACE(6,"Momentum::Error()");

    // Solve momentum equation only for interior walls
    assert(v.i>0);
    
    vd error=zeros();
    const vd& rhoti=v.rho().tdata();

    error+=Wddt*DDTfd*v.mL()();

    // Pressure right
    error+=Wpi*v.p()();
    d vSf=v.vSf;
    d vSfL=v.left()->vSf;
    // Pressure left
    error+=Wpim1*v.left()->p()();

    error+=v.mu()();
    error-=v.left()->mu()();
    // (m)^2 in time domain at i

    // Drag term
    #ifndef NODRAG
    assert(t);
    error+=Wddt*t->getDragResistance().drag(v);
    #endif
    // (Boundary) source term
    error+=v.msource();
    return error;
  }
  JacRow Momentum::jac() const {
    TRACE(6,"Momentum::jac()");
    JacRow jac(dofnr,9);

    // Solve momentum equation only for interior walls
    assert(v.i>0);

    // Time-derivative of mass flow
    jac+=JacCol(v.mL(),Wddt*DDTfd);

    jac+=JacCol(v.p(),Wpi*eye());
    jac+=JacCol(v.left()->p(),Wpim1*eye());
    jac+=JacCol(v.mu(),eye());
    jac+=JacCol(v.left()->mu(),-eye());

 
    #ifndef NODRAG
    jac+=JacCol(v.mL(),Wddt*t->getDragResistance().dm(v));
    #endif
   
    return jac;
  }
  void Momentum::domg(vd & domg_) const {
    TRACE(0,"Momentum::domg()");
    // Possibly later adding drag->domg();
    const us& Ns=v.gc->Ns();
    vd domg_full;//=v.gc->DDTfd*Wddt*fDFT*(v.rho().tdata()%v.U().tdata())/v.gc->getomg();
    // domg_.subvec(dofnr+1,dofnr+2)=domg_full.subvec(1,2);
    domg_.subvec(dofnr,dofnr+v.gc->Ns()-1)=domg_full;
    TRACE(0,"Momentum::domg() done");
  }

  #define eye (arma::eye(v.gc->Ns(),v.gc->Ns()))
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
  } // namespace LEFT
  namespace RIGHT {
   static std::tuple<d,d> weightfactors(const Cell& v){
      d xR=v.xR;
      d vxm1=v.vx;
      d vxm2=v.left()->vx;

      // Compute weight factors
      d wRNm1=(vxm2-xR)/(vxm2-vxm1);
      d wRNm2=(xR-vxm1)/(vxm2-vxm1);

      VARTRACE(25,wRNm1);
      VARTRACE(25,wRNm2);
      return std::make_tuple(wRNm1,wRNm2);
    }
   // Extrapolate momentum flow left side  
   static vd extrapolateMomentumFlow(const Cell& v){
      d wRNm1,wRNm2; std::tie(wRNm1,wRNm2)=weightfactors(v);
      return wRNm1*v.mu()()+wRNm2*v.left()->mu()();
    }
   static JacRow dExtrapolateMomentumFlow(const Cell& v){
      d wRNm1,wRNm2; std::tie(wRNm1,wRNm2)=weightfactors(v);

      JacRow jacrow(2);
      jacrow+=JacCol(v.mu(),wRNm1*eye);
      jacrow+=JacCol(v.left()->mu(),wRNm2*eye);
      return jacrow;
    }
  } // namespace RIGHT
  vd Momentum::extrapolateMomentumFlow(const Cell& v){
    TRACE(15,"Momentum::extrapolateMomentumFlow(const Cell& v)");
    assert((!v.left() && v.right()) || (!v.right() && v.left()));
    if(!v.left() && v.right())
      return LEFT::extrapolateMomentumFlow(v);
    else
      return RIGHT::extrapolateMomentumFlow(v);
  }
  JacRow Momentum::dExtrapolateMomentumFlow(const Cell& v){
    TRACE(15,"Momentum::dExtrapolateMomentumFlow(const Cell& v)");
    assert((!v.left() && v.right()) || (!v.right() && v.left()));
    if(!v.left() && v.right())
      return LEFT::dExtrapolateMomentumFlow(v);
    else
      return RIGHT::dExtrapolateMomentumFlow(v);
  }




} // namespace tube
