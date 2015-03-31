// #define TRACERPLUS 20

#include "tube.h"
#include "cell.h"
#include "weightfactors.h"
#include "jacobian.h"
#include "var.h"
#include "momentum.h"

#ifndef NODRAG
#include "drag.h"
#endif

#ifdef NODRAG
#warning Drag is turned off!
#endif

namespace tube{
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;
  using variable::var;


  void Momentum::show() const{
    cout << "----------------- Momentum equation\n";
    // cout << "Wddt   : " << Wddt << "\n";
    // cout << "WuL  : " << WuL << "\n";
    // cout << "Wu    : " << Wu << "\n";
    // cout << "WuR  : " << WuR << "\n";
    // cout << "WpL    : " << WpL << "\n";
    // cout << "WpR    : " << WpR << "\n";    
  }
  void Momentum::init()
  {
    TRACE(5,"Momentum::init(tube)");
    const Tube& t=v.getTube();
    const WeightFactors& w=v.weightFactors();
    drag=&t.getDragResistance();

    if(v.left()){
      Wddt=v.vx-v.left()->vx;;
      Wpi=w.SfL;
      Wpim1=v.SfL;
    }
  }
  vd Momentum::error() const {		// Error in momentum equation

    TRACE(6,"Momentum::Error()");
    vd error(v.gc->Ns(),fillwith::zeros);
    const vd& rhoti=v.rho().tdata();

    error+=Wddt*v.gc->DDTfd*v.rhoUL()();
    
    // Pressure terms    
    error+=Wpi*v.p()();
    if(v.left()){
      error+=Wpim1*v.left()->p()();
      const var& rhoim1=v.left()->rho();
      d vSfL=v.left()->weightFactors().vSf;

      // (rhoU)^2 in time domain at i-1
      vd rhoUsqim1t=pow(0.5*(v.left()->rhoUL().tdata()+v.rhoUL().tdata()),2);
      vd rhoUim1=fDFT*(rhoUsqim1t/(rhoim1.tdata()*vSfL));
    }

    error+=Wpi*v.p()();

    d vSf=v.weightFactors().vSf;

    // (rhoU)^2 in time domain at i
    vd rhoUsqit=pow(0.5*(v.rhoUL().tdata()+v.rhoUR().tdata()),2);
    vd rhoUi=fDFT*(rhoUsqit/(v.rho().tdata()*vSf));

    // Drag term
    assert(drag!=nullptr);
    #ifndef NODRAG
    error+=Wddt*drag->drag(v);
    #endif
    // (Boundary) source term
    error+=v.msource();
    return error;
  }
  JacRow Momentum::jac() const {
    TRACE(6,"Momentum::jac()");
    JacRow jac(dofnr,9);
    TRACE(0,"Momentum, dofnr jac:"<< dofnr);
    // jac+=drho();
    // jac+=dU();
    // jac+=dpL();
    // jac+=dpR();
    // jac+=drhoL();
    // jac+=dUL();
    // jac+=drhoR();
    // jac+=dUR();
 
    #ifndef NODRAG
    jac+=JacCol(v.U(),Wddt*drag->dUi(v));
    #endif
   
    return jac;
  }
  void Momentum::domg(vd & domg_) const {
    TRACE(0,"Momentum::domg()");
    // Possibly later adding drag->domg();
    const dmat& DDTfd=v.gc->DDTfd;
    const us& Ns=v.gc->Ns();
    vd domg_full;//=v.gc->DDTfd*Wddt*fDFT*(v.rho().tdata()%v.U().tdata())/v.gc->getomg();
    // domg_.subvec(dofnr+1,dofnr+2)=domg_full.subvec(1,2);
    domg_.subvec(dofnr,dofnr+v.gc->Ns()-1)=domg_full;
    TRACE(0,"Momentum::domg() done");
  }

  // JacCol Momentum::dU() const {
  //   TRACE(0,"Momentum::dU()");
  //   JacCol dU(v.U());
  //   dU+=Wddt*v.gc->DDTfd*fDFT*v.rho().diagt()*iDFT; // Time-derivative term
  //   dU+=2.0*Wu*fDFT*(v.rho().diagt()*v.U().diagt())*iDFT;
  //   assert(drag!=nullptr);

  //   #ifndef NODRAG
  //   dU+=Wddt*drag->dUi(v);
  //   #endif

  //   return dU;
  // }
  // JacCol Momentum::drho() const {
  //   TRACE(0,"Momentum::drho()");
  //   JacCol drho(v.rho());
  //   drho+=Wddt*v.gc->DDTfd*fDFT*v.U().diagt()*iDFT;
  //   drho+=Wu*fDFT*v.U().diagt()*v.U().diagt()*iDFT;
  //   return drho;
  // }
  // JacCol Momentum::drhoL() const {
  //   TRACE(0,"Momentum::drhoL()");
  //   JacCol drhoL(v.rhoL());
  //   drhoL+=WuL*fDFT*v.UL().diagt()*v.UL().diagt()*iDFT;
  //   return drhoL;
  // }
  // JacCol Momentum::dUL() const {
  //   TRACE(0,"Momentum::dUL()");    // Todo: add this term!;

  //   JacCol dUL(v.UL());
  //   dUL+=2.0*WuL*fDFT*v.rhoL().diagt()*v.UL().diagt()*iDFT;
  //   return dUL;
  // }
  // JacCol Momentum::drhoR() const {
  //   TRACE(0,"Momentum::dhoR()");    // Todo: add this term!;
  //   JacCol drhoR(v.rhoR());
  //   drhoR+=WuR*fDFT*v.UR().diagt()*v.UR().diagt()*iDFT;
  //   return drhoR;
  // }
  // JacCol Momentum::dUR() const {
  //   TRACE(0,"Momentum::dUR()"); // Todo: add this term!;
  //   JacCol dUR(v.UR());
  //   dUR+=2.0*WuR*fDFT*v.rhoR().diagt()*v.UR().diagt()*iDFT;
  //   return dUR;
  // }
  // JacCol Momentum::dpR() const {
  //   TRACE(0,"Momentum::dpR()");
  //   dmat I(v.gc->Ns(),v.gc->Ns(),fillwith::eye);
  //   JacCol dpR(v.pR(),WpR*I);
  //   return dpR;
  // }
  // JacCol Momentum::dpL() const {
  //   TRACE(0,"Momentum::dpR()");
  //   dmat I(v.gc->Ns(),v.gc->Ns(),fillwith::eye);
  //   JacCol dpL(v.pL(),WpL*I);
  //   return dpL;
  // }

} // namespace tube
