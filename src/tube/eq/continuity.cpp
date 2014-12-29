#include "tubevertex.h"
#include "weightfactors.h"
#include "jacobian.h"
#include "var.h"
#include "continuity.h"

#define iDFT (v.gc->iDFT)
#define fDFT (v.gc->fDFT)

namespace tube{

  using tasystem::JacCol;
  using tasystem::JacRow;
  using variable::var;
  
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
  vd Continuity::massFlow() const{
    TRACE(15,"Continuity::massFlow()");
    return fDFT*(v.rho().tdata()%v.U().tdata());
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
    error+=Wi*fDFT*(v.rho().tdata()%v.U().tdata());
    const vd& rhoL=v.rhoL().tdata();
    const vd& UL=v.UL().tdata();
    error+=WL*fDFT*(rhoL%UL);
    const vd& rhoR=v.rhoR().tdata();
    const vd& UR=v.UR().tdata();
    error+=WR*fDFT*(rhoR%UR);

    // (Boundary) source term
    error+=v.csource();
    return error;
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
    drho+=Wi*fDFT*v.U().diagt()*iDFT;
    return drho;
  }
  JacCol Continuity::dU() const {
    TRACE(0,"Continuity::dU()");
    return JacCol(v.U(),Wi*fDFT*v.rho().diagt()*iDFT);
  }
  JacCol Continuity::dUR() const {
    TRACE(0,"Continuity::dUR()");
    return JacCol(v.UR(),WR*fDFT*v.rhoR().diagt()*iDFT);
  }
  JacCol Continuity::dUL() const {
    TRACE(0,"Continuity::dUL()");
    return JacCol(v.UL(),WL*fDFT*v.rhoL().diagt()*iDFT);
  }
  JacCol Continuity::drhoR() const {
    TRACE(0,"Continuity::drhoR()");
    return JacCol(v.rhoR(),WR*fDFT*v.UR().diagt()*iDFT);
  }
  JacCol Continuity::drhoL() const {
    TRACE(0,"Continuity::drhoL()");
    return JacCol(v.rhoL(),WL*fDFT*v.UL().diagt()*iDFT);
  }
} // Namespace tube

