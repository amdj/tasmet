// #define TRACERPLUS 20
#ifdef NODRAG
#error NODRAG already defined!
#endif
// #define NODRAG

#include "momentumeq.h"
#include "tube.h"
#include "tubevertex.h"
#include "weightfactors.h"
#include "artvisco.h"
#include "jacobian.h"

namespace tube{
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;

  void Momentum::show() const{
    cout << "----------------- Momentum equation\n";
    cout << "Wddt   : " << Wddt << "\n";
    cout << "Wuim1  : " << Wuim1 << "\n";
    cout << "Wui    : " << Wui << "\n";
    cout << "Wuip1  : " << Wuip1 << "\n";
    cout << "WpL    : " << WpL << "\n";
    cout << "WpR    : " << WpR << "\n";    
    
  }
  void Momentum::init(const WeightFactors& w,const Tube& t)
  {
    TRACE(5,"Momentum::init(tube)");
    drag=&t.getDragResistance();

    // Always the same
    WpL=-w.vSf;
    WpR= w.vSf;

    Wddt=w.vVf/w.vSf;

    if(v.left()&&v.right());
    // This one should be correct
    Wuim1=-w.wLl/w.SfL;
    Wui=(w.wRl/w.SfR-w.wLr/w.SfL);
    Wuip1=w.wRr/w.SfR;

    // But this is also possible

    // m.Wuim1=-wLl/vSfL;
    // m.Wui=(wRl/w.vSf-wLr/w.vSf);
    // m.Wuip1=wRr/vSfR;





  }
  
  vd Momentum::error() const {		// Error in momentum equation

    TRACE(6,"Momentum::Error()");
    vd error(v.gc->Ns(),fillwith::zeros);
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;
    const vd& rhoti=v.rho().tdata();
    const vd& Uti=v.U().tdata();
    error+=Wddt*v.gc->DDTfd*v.gc->fDFT*(Uti%rhoti);
    error+=Wui*v.gc->fDFT*(rhoti%Uti%Uti);
    
    // Pressure terms    
    error+=WpL*v.pL()();
    error+=WpR*v.pR()();    

    if(v.left()!=NULL){
      const vd& rhotim1=v.rhoL().tdata();
      const vd& Utim1=v.UL().tdata();
      error+=Wuim1*v.gc->fDFT*(rhotim1%Utim1%Utim1);
      // Pressure term
    }
    if(v.right()!=NULL){
      const vd& Utip1=v.UR().tdata();
      const vd& rhotip1=v.rhoR().tdata();
      error+=Wuip1*v.gc->fDFT*(rhotip1%Utip1%Utip1);
    }

    // Drag term
    assert(drag!=NULL);
    #ifndef NODRAG
    error+=Wddt*drag->drag(v);
    #else
    if(v.i==0)
      TRACE(25,"Drag is turned off!");
    #endif
    // (Boundary) source term
    error+=v.msource();
    return error;
  }
  JacRow Momentum::jac() const {
    TRACE(6,"Momentum::jac()");
    JacRow jac(dofnr,9);
    TRACE(0,"Momentum, dofnr jac:"<< dofnr);
    jac+=drhoi();
    jac+=dUi();
    jac+=dpL();
    jac+=dpR();
    if(v.left()){
      jac+=drhoim1();
      jac+=dUim1();
    }
    if(v.right()){
      jac+=drhoip1();
      jac+=dUip1();
    }
    return jac;
  }
  void Momentum::domg(vd & domg_) const {
    TRACE(0,"Momentum::domg()");
    // Possibly later adding drag->domg();
    const dmat& DDTfd=v.gc->DDTfd;
    const dmat& fDFT=v.gc->fDFT;
    const us& Ns=v.gc->Ns();
    vd domg_full=v.gc->DDTfd*Wddt*fDFT*(v.rho().tdata()%v.U().tdata())/v.gc->getomg();
    // domg_.subvec(dofnr+1,dofnr+2)=domg_full.subvec(1,2);
    domg_.subvec(dofnr,dofnr+v.gc->Ns()-1)=domg_full;
    TRACE(0,"Momentum::domg() done");
  }
  JacCol Momentum::dUi() const {
    TRACE(0,"Momentum::dUi()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;

    JacCol dUi(v.U());
    // dUi+=vVf*tube.drag.dUi(i)/vSf;		       // Drag term
    dUi+=Wddt*v.gc->DDTfd*v.gc->fDFT*v.rho().diagt()*v.gc->iDFT; // Time-derivative term
    dUi+=2.0*Wui*v.gc->fDFT*(v.rho().diagt()*v.U().diagt())*v.gc->iDFT;
    // Artificial viscosity terms
    return dUi;
  }
  JacCol Momentum::drhoi() const {
    TRACE(0,"Momentum::drhoi()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;
    JacCol drhoi(v.rho());
    drhoi+=Wddt*v.gc->DDTfd*v.gc->fDFT*v.U().diagt()*v.gc->iDFT;
    drhoi+=Wui*v.gc->fDFT*v.U().diagt()*v.U().diagt()*v.gc->iDFT;
    return drhoi;
  }
  JacCol Momentum::drhoim1() const {
    TRACE(0,"Momentum::drhoim1()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;
    JacCol drhoim1(v.rhoL());
    drhoim1+=Wuim1*v.gc->fDFT*v.UL().diagt()*v.UL().diagt()*v.gc->iDFT;
    return drhoim1;
  }
  JacCol Momentum::dUim1() const {
    TRACE(0,"Momentum::dUim1()");    // Todo: add this term!;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;

    JacCol dUim1(v.UL());
    dUim1+=2.0*Wuim1*v.gc->fDFT*v.rhoL().diagt()*v.UL().diagt()*v.gc->iDFT;
    // Artificial viscosity terms
    
    return dUim1;
  }
  JacCol Momentum::drhoip1() const {
    TRACE(0,"Momentum::dhoip1()");    // Todo: add this term!;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;
    JacCol drhoip1(v.rhoR());
    drhoip1+=Wuip1*v.gc->fDFT*v.UR().diagt()*v.UR().diagt()*v.gc->iDFT;
    return drhoip1;
  }
  JacCol Momentum::dUip1() const {
    TRACE(0,"Momentum::dUip1()"); // Todo: add this term!;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;

    JacCol dUip1(v.UR());
    dUip1+=2.0*Wuip1*v.gc->fDFT*v.rhoR().diagt()*v.UR().diagt()*v.gc->iDFT;
    return dUip1;
  }
  JacCol Momentum::dpR() const {
    TRACE(0,"Momentum::dpR()");

    dmat I(v.gc->Ns(),v.gc->Ns(),fillwith::eye);
    JacCol dpR(v.pR(),WpR*I);
    return dpR;
  }
  JacCol Momentum::dpL() const {
    TRACE(0,"Momentum::dpR()");

    dmat I(v.gc->Ns(),v.gc->Ns(),fillwith::eye);
    JacCol dpL(v.pL(),WpL*I);
    return dpL;
  }
  // JacCol Momentum::dUip2() const {
  //   TRACE(0,"Momentum:dUip2()");
  //   // TRACE(50,"i:"<<i);
  //   dmat dUip2=v.zero;
  //   #ifdef MOM_VISCOSITY
  //   // if(v.left==NULL)
  //     // dUip2+=-v.gc->rho0*d_l();
  //   #endif
  //   return dUip2;
  // }
  // JacCol Momentum::dUim2() const {
  //   dmat dUim2=v.zero;
  //   #ifdef MOM_VISCOSITY
  //   // if(v.right()==NULL )
  //     // dUim2+=-v.gc->rho0*d_r();
  //   #endif
  //   return dUim2;
  // }  

} // namespace tube
