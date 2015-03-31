#include "cell.h"
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
  }
  void Continuity::init(){
    TRACE(8,"Continuity::init(tube)");
    Wddt=v.vVf;
  } // init
  JacRow Continuity::jac() const{
    TRACE(6,"Continuity::jac()");
    JacRow jac(dofnr,3);
    // TRACE(0,"Continuity, dofnr jac:"<< dofnr);
    jac.addCol(drho());    
    jac.addCol(dmL());
    jac.addCol(dmR());
    return jac;
  }
  vd Continuity::error() const {	
    TRACE(6,"Continuity::Error()");
    vd error(v.gc->Ns(),fillwith::zeros);
    error+=Wddt*v.gc->DDTfd*v.rho()();
    error+=v.mR();
    error-=v.mL();

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
  JacCol Continuity::dmR() const {
    TRACE(0,"Continuity::drhoR()");
    return JacCol(v.mR(),eye());
  }
  JacCol Continuity::dmL() const {
    TRACE(0,"Continuity::drhoL()");
    return JacCol(v.mL(),eye());
  }
} // Namespace tube

