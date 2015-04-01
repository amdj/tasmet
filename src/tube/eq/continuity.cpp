#include "cell.h"
#include "jacobian.h"
#include "var.h"
#include "continuity.h"

#define iDFT (v.gc->iDFT)
#define fDFT (v.gc->fDFT)
#define DDTfd (v.gc->DDTfd)

namespace tube{

  using tasystem::JacCol;
  using tasystem::JacRow;
  using variable::var;
  
  void Continuity::show() const {
    cout << "-------------- Continuity equation\n";
    cout << "Wddt : " << Wddt << "\n";
  }
  Continuity::~Continuity(){TRACE(0,"~Continuity()");}
  void Continuity::init(){
    TRACE(8,"Continuity::init(tube)");
    Wddt=v.vVf;
  } // init
  JacRow Continuity::jac() const{
    TRACE(6,"Continuity::jac()");
    JacRow jac(dofnr,3);
    // TRACE(0,"Continuity, dofnr jac:"<< dofnr);
    jac+=JacCol(v.rho(),Wddt*DDTfd);    
    jac+=JacCol(v.mR(),eye());
    jac+=JacCol(v.mL(),-eye());
    return jac;
  }
  vd Continuity::error() const {	
    TRACE(6,"Continuity::Error()");
    vd error(v.gc->Ns(),fillwith::zeros);
    error+=Wddt*DDTfd*v.rho()();
    error+=v.mR()();
    error-=v.mL()();

    // (Boundary) source term
    error+=v.csource();
    return error;
  }
  void Continuity::domg(vd& domg_) const{
    TRACE(0,"Continuity::domg()");
    vd domg_full=Wddt*DDTfd*v.rho()()/v.gc->getomg();
    // domg_.subvec(dofnr+1,dofnr+2)=domg_full.subvec(1,2); 
    domg_.subvec(dofnr,dofnr+v.gc->Ns()-1)=domg_full;     
  }
} // Namespace tube

