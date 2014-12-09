// #define TRACERPLUS 20
#ifdef NODRAG
#error NODRAG already defined!
#endif
// #define NODRAG

#include "momentumeq.h"
#include "tube.h"
#include "tubevertex.h"
#include "weightfactors.h"
#include "jacobian.h"
#include "var.h"

namespace tube{
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;
  using variable::var;
  void Momentum::show() const{
    cout << "----------------- Momentum equation\n";
    cout << "Wddt   : " << Wddt << "\n";
    cout << "WuL  : " << WuL << "\n";
    cout << "Wu    : " << Wu << "\n";
    cout << "WuR  : " << WuR << "\n";
    cout << "WpL    : " << WpL << "\n";
    cout << "WpR    : " << WpR << "\n";    
    
  }
  void Momentum::init()
  {
    TRACE(5,"Momentum::init(tube)");
    const Tube& t=v.getTube();
    const WeightFactors& w=w.weightFactors();
    drag=t.getDragResistance();

    // Always the same
    WpL=-w.vSf;
    WpR= w.vSf;

    Wddt=w.vVf/w.vSf;

    if(v.left()){
      WuL+=-w.wLl/w.vSfL;
      Wu+=-w.wLr/w.vSf;
    }
    else{
      WuL+=-1/w.SfL;
    }
    if(v.right()){
      WuR+=w.wRr/w.vSfR;
      Wu+=w.wRl/w.vSf;
    }
    else{
      WuR+=1/w.SfR;
    }      
    // For paper 2, we have used:

    // WuL=-w.wLl/w.vSfL;
    // Wu=(w.wRl/w.SfR-w.wLr/w.SfL);
    // WuR=w.wRr/w.SfR;
  }
  
  vd Momentum::error() const {		// Error in momentum equation

    TRACE(6,"Momentum::Error()");
    vd error(v.gc->Ns(),fillwith::zeros);
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;
    const vd& rhoti=v.rho().tdata();
    const vd& Uti=v.U().tdata();
    error+=Wddt*v.gc->DDTfd*v.gc->fDFT*(Uti%rhoti);
    error+=Wu*v.gc->fDFT*(rhoti%Uti%Uti);
    
    // Pressure terms    
    error+=WpL*v.pL()();
    error+=WpR*v.pR()();    

    if(v.left()){
      const vd& rhotL=v.rhoL().tdata();
      const vd& UtL=v.UL().tdata();
      error+=WuL*v.gc->fDFT*(rhotL%UtL%UtL);
      // Pressure term
    }
    if(v.right()){
      const vd& UtR=v.UR().tdata();
      const vd& rhotR=v.rhoR().tdata();
      error+=WuR*v.gc->fDFT*(rhotR%UtR%UtR);
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
    bjac+=drho();
    jac+=dUi();
    jac+=dpL();
    jac+=dpR();
    jac+=drhoL();
    jac+=dUL();
    jac+=drhoR();
    jac+=dUR();
    
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
  var Momentum::momentumFlow() const {
    TRACE(3,"Momentum::momentumFlow()");
    const weightFactors& w=v.weightFactors();
    vd Ut=v.U().tdata();
    vd rhot=v.rho().tdata();
    var res(v.gc);
    res.settdata(v.gc->fDFT*(rhot%Ut%Ut/w.vSf));
    return res;
  }
  vd Momentum::extrapolateMomentumFlow() const{
    vd momentumFlow;
    const WeightFactors& w=v.weightFactors();
    if(!v.left()){
      const vd& rhoR=v.rhoR().tdata();
      const vd& UR=v.UR().tdata();
      momentumFlow=w.wL1*v.gc->fDFT*(rhoR%UR%UR/w.vSfR);
      momentumFlow+=w.wL0*momentumFlow()();
      return momentumFlow;
    }
    if(!v.right()){
      const vd& rhoL=v.rhoL().tdata();
      const vd& UL=v.UL().tdata();
      momentumFlow=w.wRNm2*v.gc->fDFT*(rhoL%UL%UL/w.vSfL);
      momentumFlow+=w.wRNm1*momentumFlow()();
      return momentumFlow;
    }
    else{
      WARN("SOMETHING REALLY WRONG! Momentumflow tried to be extrapolated on a non-boundary vertex");
    }
  }
  JacRow Momentum::dExtrapolateMassFlow() const{
    const weightFactors& w=v.weightFactors();
    JacRow jacrow(-1,4);
    dmat Utd=v.U().diagt();
    dmat rhotd=v.rho().diagt();

    if(!v.left()){
      JacCol dU(v.U(),2.0*w.wL0/w.vSf*v.gc->fDFT\
                *rhotd*Utd*gc->iDFT);
      JacCol drho(v.rho(),(w.wL0/w.vSf)*v.gc->fDFT\
                  *(Utd%Utd)*gc->iDFT);

      dmat UtdR=v.UR().diagt();
      dmat rhotdR=v.rhoR().diagt();
      JacCol dUR(v.UR(),2.0*(w.wL1/w.vSfR)v.gc->fDFT\
                 *rhotdR*UtdR*gc->iDFT);
      JacCol drhoR(v.rhoR(),(w.wL1/w.vSfR)*v.gc->fDFT\
                   *UtdR*UtdR*gc->iDFT);
      jacrow+=dU+=drho+=dUR+=drhoR;
      return jacrow;
    }
    if(!v.right()){
      JacCol dU(v.U(),2.0*w.wRNm1/w.vSf*v.gc->fDFT\
                *rhotd*Utd*gc->iDFT);
      JacCol drho(v.rho(),(w.wRNm1/w.vSf)*v.gc->fDFT\
                  *(Utd%Utd)*gc->iDFT);

      dmat UtdL=v.UL().diagt();
      dmat rhotdL=v.rhoL().diagt();
      JacCol dUL(v.UL(),2.0*(w.wRNm2/w.vSfL)v.gc->fDFT\
                 *rhotdL*UtdL*gc->iDFT);
      JacCol drhoL(v.rhoL(),(w.wLNm2/w.vSfL)*v.gc->fDFT\
                   *UtdL*UtdL*gc->iDFT);
      jacrow+=dU+=drho+=dUL+=drhoL;
      return jacrow;
    }
  }

  JacCol Momentum::dU() const {
    TRACE(0,"Momentum::dU()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;

    JacCol dU(v.U());
    // dUi+=vVf*tube.drag.dUi(i)/vSf;		       // Drag term
    dU+=Wddt*v.gc->DDTfd*v.gc->fDFT*v.rho().diagt()*v.gc->iDFT; // Time-derivative term
    dU+=2.0*Wu*v.gc->fDFT*(v.rho().diagt()*v.U().diagt())*v.gc->iDFT;
    return dU;
  }
  JacCol Momentum::drho() const {
    TRACE(0,"Momentum::drho()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;
    JacCol drho(v.rho());
    drho+=Wddt*v.gc->DDTfd*v.gc->fDFT*v.U().diagt()*v.gc->iDFT;
    drho+=Wu*v.gc->fDFT*v.U().diagt()*v.U().diagt()*v.gc->iDFT;
    return drho;
  }
  JacCol Momentum::drhoL() const {
    TRACE(0,"Momentum::drhoL()");
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;
    JacCol drhoL(v.rhoL());
    drhoL+=WuL*v.gc->fDFT*v.UL().diagt()*v.UL().diagt()*v.gc->iDFT;
    return drhoL;
  }
  JacCol Momentum::dUL() const {
    TRACE(0,"Momentum::dUL()");    // Todo: add this term!;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;

    JacCol dUL(v.UL());
    dUL+=2.0*WuL*v.gc->fDFT*v.rhoL().diagt()*v.UL().diagt()*v.gc->iDFT;
    return dUL;
  }
  JacCol Momentum::drhoR() const {
    TRACE(0,"Momentum::dhoR()");    // Todo: add this term!;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;
    JacCol drhoR(v.rhoR());
    drhoR+=WuR*v.gc->fDFT*v.UR().diagt()*v.UR().diagt()*v.gc->iDFT;
    return drhoR;
  }
  JacCol Momentum::dUR() const {
    TRACE(0,"Momentum::dUR()"); // Todo: add this term!;
    const dmat& fDFT=v.gc->fDFT;
    const dmat& iDFT=v.gc->iDFT;

    JacCol dUR(v.UR());
    dUR+=2.0*WuR*v.gc->fDFT*v.rhoR().diagt()*v.UR().diagt()*v.gc->iDFT;
    return dUR;
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

} // namespace tube
