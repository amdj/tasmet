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
    #ifndef NODRAG
    drag=&t.getDragResistance();
    #endif
    if(v.left()){
      Wddt=v.vx-v.left()->vx;;
      Wpi=v.SfL;
      Wpim1=-v.SfL;
      Wmim1=-1/(v.left()->vSf);
      Wmi=1/(v.left()->vSf);
    }
  }
  vd Momentum::error() const {		// Error in momentum equation
    TRACE(6,"Momentum::Error()");

    // Solve momentum equation only for interior walls
    assert(v.i>0);

    
    vd error=zeros();
    const vd& rhoti=v.rho().tdata();

    error+=Wddt*DDTfd*v.mL()();

    // Pressure terms    
    error+=Wpi*v.p()();
    d vSf=v.vSf;

    // (m)^2 in time domain at i
    vd m_sq_i_td=pow(0.5*(v.mL().tdata()+v.mR().tdata()),2);
    vd mi=fDFT*(Wmi*m_sq_i_td/(v.rho().tdata()));

    // Minus what comes in
    error+=Wpim1*v.left()->p()();
    const var& rhoim1=v.left()->rho();
    d vSfL=v.left()->vSf;

    // (m)^2 in time domain at i-1
    vd m_sq_im1_td=pow(0.5*(v.left()->mL().tdata()+v.mL().tdata()),2);
    vd mim1=fDFT*(Wmim1*m_sq_im1_td/(rhoim1.tdata()));

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

    // Solve momentum equation only for interior walls
    assert(v.i>0);


    // Time-derivative of mass flow
    jac+=JacCol(v.mL(),Wddt*DDTfd);

    jac+=JacCol(v.p(),Wpi*eye());
    jac+=JacCol(v.left()->p(),Wpim1*eye());

    // Mass flow at node i-1 in time domain [kg/s]
    vd m_im1_td=0.5*(v.left()->mL().tdata()+v.mL().tdata());
    // Mass flow at node i-1 in time domain squared
    vd m_sq_im1_td=pow(m_im1_td,2);

    // Same for at node i
    vd m_i_td=0.5*(v.mL().tdata()+v.mR().tdata());
    vd m_sq_i_td=pow(m_i_td,2);    

    const var& rhoL=v.left()->rho();
    const var& rho=v.rho();    
    jac+=JacCol(v.rho(),-Wmi*fDFT*diagmat(m_sq_i_td/pow(rho.tdata(),2))*iDFT);
    jac+=JacCol(v.left()->rho(),-Wmim1*fDFT*diagmat(m_sq_im1_td/pow(rhoL.tdata(),2))*iDFT);                

    jac+=JacCol(v.mL(),fDFT*diagmat(Wmi*m_i_td/rho.tdata())*iDFT);
    jac+=JacCol(v.mL(),fDFT*diagmat(Wmim1*m_im1_td/rho.tdata())*iDFT);

    jac+=JacCol(v.left()->mL(),fDFT*diagmat(Wmim1*m_im1_td/rho.tdata())*iDFT);    
    jac+=JacCol(v.mR(),fDFT*diagmat(Wmi*m_i_td/rho.tdata())*iDFT);    
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
  // JacCol Momentum::dpi() const {
  //   TRACE(0,"Momentum::dpi()");
  //   return JacCol(v.p(),Wpi*eye());
  // }
  // JacCol Momentum::dpim1() const {
  //   TRACE(0,"Momentum::dpim1()");
  //   return JacCol(v.left()->p(),Wpim1*eye());
  // }
  // JacCol Momentum::dpL() const {
  //   TRACE(0,"Momentum::dpR()");
  //   dmat I(v.gc->Ns(),v.gc->Ns(),fillwith::eye);
  //   JacCol dpL(v.pL(),WpL*I);
  //   return dpL;
  // }

} // namespace tube
