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
    const Tube& t=v.getTube();
    #ifndef NODRAG
    drag=&t.getDragResistance();
    #endif
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
    assert(drag!=nullptr);
    error+=Wddt*drag->drag(v);
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
    jac+=JacCol(v.mL(),Wddt*drag->dm(v));
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

} // namespace tube
