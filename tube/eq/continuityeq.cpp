#include "continuityeq.h"
#include "tubevertex.h"
#include "weightfactors.h"
#include "jacobian.h"

namespace tube{

  using tasystem::JacCol;
  using tasystem::JacRow;
  
  void Continuity::show() const {
    cout << "-------------- Continuity equation\n";
    cout << "Wddt : " << Wddt << "\n";
    cout << "WL : " << WL << "\n";
    cout << "Wi   : " << Wi  << "\n";
    cout << "WR : " << WR << "\n";    
  }
  void Continuity::init(){
    TRACE(8,"Continuity::init(tube)");
    const WeightFactors& w=v.weightFactors();
    Wddt=w.vVf;

    if(v.left()&&v.right()){
      WL=-w.wLl;
      Wi=w.wRl-w.wLr;
      WR=w.wRr;
    }
    else if(v.right()){         // Leftmost vertex
      // Leftmost vertex
      Wi=w.wRl;
      WR=w.wRr;
      WL=-1;
    }
    else if(v.left()){          // Rightmost vertex
      // Leftmost vertex
      Wi=-w.wLr;
      WR=1;
      WL=-w.wLl;
    }
  } // init
  var Continuity::massFlow() const{
    return v.rho()*v.U();
  }
  JacRow Continuity::jac() const{
    TRACE(6,"Continuity::jac()");
    JacRow jac(dofnr,6);
    TRACE(0,"Continuity, dofnr jac:"<< dofnr);
    jac.addCol(drho());
    jac.addCol(dU());
    jac.addCol(drhoL());    
    jac.addCol(dUL());
    jac.addCol(drhoR());
    jac.addCol(dUR());

    return jac;
  }
  vd Continuity::error() const {	
    TRACE(6,"Continuity::Error()");
    vd error(v.gc->Ns(),fillwith::zeros);
    error+=Wddt*v.gc->DDTfd*v.rho()();
    error+=Wi*v.gc->fDFT*(v.rho().tdata()%v.U().tdata());
    if(v.left()){
      const vd& rhoL=v.rhoL().tdata();
      const vd& UL=v.UL().tdata();
      error+=WL*v.gc->fDFT*(rhoL%UL);
    }
    // TRACE(10,"Right:,"<<v.right());    
    if(v.right()){
      // Standard implementation of a no-slip (wall) boundary
      // condition
      const vd& rhoR=v.rhoR().tdata();
      const vd& UR=v.UR().tdata();
      error+=WR*v.gc->fDFT*(rhoR%UR);
    }
    // (Boundary) source term
    error+=v.csource();
    return error;
  }
  vd Continuity::extrapolateMassFlow() const{
    vd massFlow;
    const WeightFactors& w=v.weightFactors();
    if(!v.left()){
      const vd& rhoR=v.rhoR().tdata();
      const vd& UR=v.UR().tdata();
      massFlow=w.wL1*v.gc->fDFT*(rhoR%UR);
      massFlow+=w.wL0*massFlow();
      return massFlow;
    }
    if(!v.right()){
      const vd& rhoL=v.rhoL().tdata();
      const vd& UL=v.UL().tdata();
      massFlow=w.wRNm2*v.gc->fDFT*(rhoL%UL);
      massFlow+=w.wRNm1*massFlow();
      return massFlow;
    }
    else{
      WARN("SOMETHING REALLY WRONG! Massflow tried to be extrapolated on a non-boundary vertex");
    }
  }
  JacRow Continuity::dExtrapolateMassFlow() const{
    const weightFactors& w=v.weightFactors();
    JacRow jacrow(-1,4);
    if(!v.left()){
      JacCol dU(v.U(),w.wL0*v.gc->fDFT*v.rho().diagt()*gc->iDFT);
      JacCol drho(v.rho(),w.wL0*v.gc->fDFT*v.U().diagt()*gc->iDFT);
      JacCol dUR(v.UR(),w.wL1*v.gc->fDFT*v.rhoR().diagt()*gc->iDFT);
      JacCol drhoR(v.rhoR(),w.wL1*v.gc->fDFT*v.UR().diagt()*gc->iDFT);
      jacrow+=dU+=drho+=dUR+=drhoR;
      return jacrow;
    }
    if(!v.right()){
      JacCol dU(v.U(),w.wRNm1*v.gc->fDFT*v.rho().diagt()*gc->iDFT);
      JacCol drho(v.rho(),w.wRNm1*v.gc->fDFT*v.U().diagt()*gc->iDFT);
      JacCol dUL(v.UL(),w.wRNm2*v.gc->fDFT*v.rhoL().diagt()*gc->iDFT);
      JacCol drhoL(v.rhoL(),w.wRNm2*v.gc->fDFT*v.UL().diagt()*gc->iDFT);
      jacrow+=dU+=drho+=dUL+=drhoL;
      return jacrow;
    }
  }
  void Continuity::domg(vd& domg_) const{
    TRACE(0,"Continuity::domg()");
    vd domg_full=Wddt*v.gc->DDTfd*v.rho()()/v.gc->getomg();
    // domg_.subvec(dofnr+1,dofnr+2)=domg_full.subvec(1,2); 
    domg_.subvec(dofnr,dofnr+v.gc->Ns()-1)=domg_full;     
  }
  JacCol Continuity::drho() const {
    TRACE(0,"Continuity::drho()");
    JacCol drho(v.rho(),Wddt*v.gc->DDTfd);		// Initialize and add first term
    drho+=Wi*v.gc->fDFT*v.U().diagt()*v.gc->iDFT;
    return drho;
  }
  JacCol Continuity::dU() const {
    TRACE(0,"Continuity::dUi()");
    return JacCol(v.U(),Wi*v.gc->fDFT*v.rho().diagt()*v.gc->iDFT);
  }
  JacCol Continuity::dUR() const {
    TRACE(0,"Continuity::dUR()");
    return JacCol(v.UR(),WR*v.gc->fDFT*v.rhoR().diagt()*v.gc->iDFT);

  }
  JacCol Continuity::dUL() const {
    TRACE(0,"Continuity::dUL()");
    return JacCol(v.UL(),WL*v.gc->fDFT*v.rhoL().diagt()*v.gc->iDFT);
  }
  JacCol Continuity::drhoR() const {
    TRACE(0,"Continuity::drhoR()");
    return JacCol(v.rhoR(),WR*v.gc->fDFT*v.UR().diagt()*v.gc->iDFT);
  }
  JacCol Continuity::drhoL() const {
    TRACE(0,"Continuity::drhoL()");
    return JacCol(v.rhoL(),WL*v.gc->fDFT*v.UL().diagt()*v.gc->iDFT);
  }
} // Namespace tube

















