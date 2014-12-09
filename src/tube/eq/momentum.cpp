// #define TRACERPLUS 20
#ifdef NODRAG
#error NODRAG already defined!
#endif
// #define NODRAG

#include "tube.h"
#include "tubevertex.h"
#include "weightfactors.h"
#include "jacobian.h"
#include "var.h"
#include "momentum.h"

#define iDFT (v.gc->iDFT)
#define fDFT (v.gc->fDFT)


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
    const WeightFactors& w=v.weightFactors();
    drag=&t.getDragResistance();

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
    const vd& rhoti=v.rho().tdata();
    const vd& Uti=v.U().tdata();
    error+=Wddt*v.gc->DDTfd*fDFT*(Uti%rhoti);
    error+=Wu*fDFT*(rhoti%Uti%Uti);
    
    // Pressure terms    
    error+=WpL*v.pL()();
    error+=WpR*v.pR()();    

    if(v.left()){
      const vd& rhotL=v.rhoL().tdata();
      const vd& UtL=v.UL().tdata();
      error+=WuL*fDFT*(rhotL%UtL%UtL);
      // Pressure term
    }
    if(v.right()){
      const vd& UtR=v.UR().tdata();
      const vd& rhotR=v.rhoR().tdata();
      error+=WuR*fDFT*(rhotR%UtR%UtR);
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
    jac+=drho();
    jac+=dU();
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
    const us& Ns=v.gc->Ns();
    vd domg_full=v.gc->DDTfd*Wddt*fDFT*(v.rho().tdata()%v.U().tdata())/v.gc->getomg();
    // domg_.subvec(dofnr+1,dofnr+2)=domg_full.subvec(1,2);
    domg_.subvec(dofnr,dofnr+v.gc->Ns()-1)=domg_full;
    TRACE(0,"Momentum::domg() done");
  }
  vd Momentum::momentumFlow() const {
    TRACE(3,"Momentum::momentumFlow()");
    const WeightFactors& w=v.weightFactors();
    vd Ut=v.U().tdata();
    vd rhot=v.rho().tdata();
    var res(v.gc);
    return fDFT*(rhot%Ut%Ut/w.vSf);
  }
  vd Momentum::extrapolateMomentumFlow() const{
    vd emomentumFlow;
    const WeightFactors& w=v.weightFactors();
    if(!v.left()){
      const vd& rhoR=v.rhoR().tdata();
      const vd& UR=v.UR().tdata();
      emomentumFlow=w.wL1*v.left()->momentum().momentumFlow();
      emomentumFlow+=w.wL0*momentumFlow();
      return emomentumFlow;
    }
    if(!v.right()){
      const vd& rhoL=v.rhoL().tdata();
      const vd& UL=v.UL().tdata();
      emomentumFlow=w.wRNm2*v.right()->momentum().momentumFlow();
      emomentumFlow+=w.wRNm1*momentumFlow();
      return emomentumFlow;
    }
    else{
      WARN("SOMETHING REALLY WRONG! Momentumflow tried to be extrapolated on a non-boundary vertex");
      abort();
      return vd(0,fillwith::zeros);
    }
  }
  JacRow Momentum::dExtrapolateMomentumFlow() const{
    const WeightFactors& w=v.weightFactors();
    JacRow jacrow(-1,4);
    dmat Utd=v.U().diagt();
    dmat rhotd=v.rho().diagt();

    if(!v.left()){
      JacCol dU(v.U(),2.0*(w.wL0/w.vSf)*fDFT  \
                *rhotd*Utd*iDFT);
      JacCol drho(v.rho(),(w.wL0/w.vSf)*fDFT\
                  *(Utd%Utd)*iDFT);

      dmat UtdR=v.UR().diagt();
      dmat rhotdR=v.rhoR().diagt();
      JacCol dUR(v.UR(),2.0*(w.wL1/w.vSfR)*fDFT\
                 *rhotdR*UtdR*iDFT);
      JacCol drhoR(v.rhoR(),(w.wL1/w.vSfR)*fDFT\
                   *UtdR*UtdR*iDFT);
      ((((jacrow+=dU)+=drho)+=dUR)+=drhoR);
    }
    if(!v.right()){
      JacCol dU(v.U(),2.0*w.wRNm1/w.vSf*fDFT\
                *rhotd*Utd*iDFT);
      JacCol drho(v.rho(),(w.wRNm1/w.vSf)*fDFT\
                  *(Utd%Utd)*iDFT);

      dmat UtdL=v.UL().diagt();
      dmat rhotdL=v.rhoL().diagt();
      JacCol dUL(v.UL(),2.0*(w.wRNm2/w.vSfL)*fDFT\
                 *rhotdL*UtdL*iDFT);
      JacCol drhoL(v.rhoL(),(w.wRNm2/w.vSfL)*fDFT\
                   *UtdL*UtdL*iDFT);
      ((((jacrow+=dU)+=drho)+=dUL)+=drhoL);
    }
    else{
      WARN("SOMETHING REALLY WRONG! dMomentumflow tried to be extrapolated on a non-boundary vertex");
      abort();
    }
    return jacrow;
  }

  JacCol Momentum::dU() const {
    TRACE(0,"Momentum::dU()");
    JacCol dU(v.U());
    dU+=Wddt*v.gc->DDTfd*fDFT*v.rho().diagt()*iDFT; // Time-derivative term
    dU+=2.0*Wu*fDFT*(v.rho().diagt()*v.U().diagt())*iDFT;
    assert(drag!=NULL);

    #ifndef NODRAG
    dU+=Wddt*drag->dUi(v);
    #endif

    return dU;
  }
  JacCol Momentum::drho() const {
    TRACE(0,"Momentum::drho()");
    JacCol drho(v.rho());
    drho+=Wddt*v.gc->DDTfd*fDFT*v.U().diagt()*iDFT;
    drho+=Wu*fDFT*v.U().diagt()*v.U().diagt()*iDFT;
    return drho;
  }
  JacCol Momentum::drhoL() const {
    TRACE(0,"Momentum::drhoL()");
    JacCol drhoL(v.rhoL());
    drhoL+=WuL*fDFT*v.UL().diagt()*v.UL().diagt()*iDFT;
    return drhoL;
  }
  JacCol Momentum::dUL() const {
    TRACE(0,"Momentum::dUL()");    // Todo: add this term!;

    JacCol dUL(v.UL());
    dUL+=2.0*WuL*fDFT*v.rhoL().diagt()*v.UL().diagt()*iDFT;
    return dUL;
  }
  JacCol Momentum::drhoR() const {
    TRACE(0,"Momentum::dhoR()");    // Todo: add this term!;
    JacCol drhoR(v.rhoR());
    drhoR+=WuR*fDFT*v.UR().diagt()*v.UR().diagt()*iDFT;
    return drhoR;
  }
  JacCol Momentum::dUR() const {
    TRACE(0,"Momentum::dUR()"); // Todo: add this term!;
    JacCol dUR(v.UR());
    dUR+=2.0*WuR*fDFT*v.rhoR().diagt()*v.UR().diagt()*iDFT;
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
