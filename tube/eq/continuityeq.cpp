#include "continuityeq.h"
#include "tubevertex.h"
#include "weightfactors.h"
#include "artvisco.h"
#include "jacobian.h"

namespace tube{

  using tasystem::JacCol;
  using tasystem::JacRow;
  
  void Continuity::show() const {
    cout << "-------------- Continuity equation\n";
    cout << "Wddt : " << Wddt << "\n";
    cout << "Wim1 : " << Wim1 << "\n";
    cout << "Wi   : " << Wi  << "\n";
    cout << "Wip1 : " << Wip1 << "\n";    
    #ifdef CONT_VISCOSITY
    cout << "Artificial viscosity turned ON for continuity equation\n";
    #else
    cout << "Artificial viscosity turned OFF for continuity equation\n";
    #endif
    
  }
  void Continuity::init(const WeightFactors& w,const Tube& t){
    TRACE(8,"Continuity::init(tube)");
    Wddt=w.vVf;

    if(v.left()&&v.right()){
      Wim1=-w.wLl;
      Wi=w.wRl-w.wLr;
      Wip1=w.wRr;
    }
    else if(v.right()){         // Leftmost vertex
      // Leftmost vertex
      Wi=w.wRl;
      Wip1=w.wRr;
      Wim1=-1;
    }
    else if(v.left()){          // Rightmost vertex
      // Leftmost vertex
      Wi=-w.wLr;
      Wip1=1;
      Wim1=-w.wLl;
    }

  } // init

  JacRow Continuity::jac() const{
    TRACE(6,"Continuity::jac()");
    JacRow jac(dofnr,6);
    TRACE(0,"Continuity, dofnr jac:"<< dofnr);
    jac.addCol(drhoi());
    jac.addCol(dUi());
    if(v.left()){
      jac.addCol(drhoim1());    
      jac.addCol(dUim1());
    }
    if(v.right()){
      jac.addCol(drhoip1());
      jac.addCol(dUip1());
    }
    return jac;
  }
  
  vd Continuity::error() const {	
    TRACE(6,"Continuity::Error()");
    vd error(v.gc->Ns(),fillwith::zeros);
    error+=Wddt*v.gc->DDTfd*v.rho()();
    error+=Wi*v.gc->fDFT*(v.rho().tdata()%v.U().tdata());
    if(v.left()){
      const vd& rhoim1=v.rhoL().tdata();
      const vd& Uim1=v.UL().tdata();
      error+=Wim1*v.gc->fDFT*(rhoim1%Uim1);
    }
    // TRACE(10,"Right:,"<<v.right());    
    if(v.right()){
      // Standard implementation of a no-slip (wall) boundary
      // condition
      const vd& rhoip1=v.rhoR().tdata();
      const vd& Uip1=v.UR().tdata();
      error+=Wip1*v.gc->fDFT*(rhoip1%Uip1);
    }
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
  JacCol Continuity::drhoi() const {
    TRACE(0,"Continuity::drhoi()");
    JacCol drhoi(v.rho(),Wddt*v.gc->DDTfd);		// Initialize and add first term
    drhoi+=Wi*v.gc->fDFT*v.U().diagt()*v.gc->iDFT;
    return drhoi;
  }
  JacCol Continuity::dUi() const {
    TRACE(0,"Continuity::dUi()");
    JacCol dUi(v.U(),Wi*v.gc->fDFT*v.rho().diagt()*v.gc->iDFT);
    return dUi;
  }
  JacCol Continuity::dUip1() const {
    TRACE(0,"Continuity::dUip1()");
    JacCol dUip1(v.UR(),Wip1*v.gc->fDFT*v.rhoR().diagt()*v.gc->iDFT);
    return dUip1;
  }
  JacCol Continuity::dUim1() const {
    TRACE(0,"Continuity::dUim1()");
    JacCol dUim1(v.UL(),Wim1*v.gc->fDFT*v.rhoL().diagt()*v.gc->iDFT);
    return dUim1;
  }
  JacCol Continuity::drhoip1() const {
    TRACE(0,"Continuity::drhoip1()");
    JacCol drhoip1(v.rhoR(),Wip1*v.gc->fDFT*v.UR().diagt()*v.gc->iDFT);
    return drhoip1;
  }
  JacCol Continuity::drhoim1() const {
    TRACE(0,"Continuity::drhoim1()");

    JacCol drhoim1(v.rhoL(),Wim1*v.gc->fDFT*v.UL().diagt()*v.gc->iDFT);
    return drhoim1;
  }
} // Namespace tube

















