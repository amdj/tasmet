#include "cell.h"
#include "jacrow.h"
#include "var.h"
#include "continuity.h"
#include "weightfactors.h"

#define iDFT (v.gc->iDFT)
#define fDFT (v.gc->fDFT)
#define DDTfd (v.gc->DDTfd)
#define Ns (v.gc->Ns())
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
    vd error(Ns,fillwith::zeros);
    error+=Wddt*DDTfd*v.rho()();
    error+=v.mR()();
    error-=v.mL()();

    // (Boundary) source term
    error+=v.csource();
    return error;
  }
  void Continuity::domg(vd& domg_) const{
    TRACE(18,"Continuity::domg()");
    vd domg_full=(Wddt*DDTfd*v.rho()())/v.gc->getomg();
    domg_.subvec(dofnr,dofnr+Ns-1)=domg_full;     
  }
  vd Continuity::extrapolateMassFlow(const Cell& v){
    TRACE(15,"Continuity::extrapolateMassFlow(const Cell& v)");
    assert((!v.left() && v.right()) || (!v.right() && v.left()));
    if(!v.left()){
      d W0,W1; std::tie(W0,W1)=BcWeightFactors(v);
      return W0*v.mR()()+W1*v.right()->mR()();
    }
    else{
      d WR1,WR2; std::tie(WR1,WR2)=BcWeightFactors(v);
      return WR1*v.mL()()+WR2*v.left()->mL()();
    }
  }
  JacRow Continuity::dExtrapolateMassFlow(const Cell& v){
    TRACE(15,"Continuity::dExtrapolateMassFlow(const Cell& v)");
    assert((!v.left() && v.right()) || (!v.right() && v.left()));
    JacRow jac(2);
    if(!v.left()){
      d W0,W1; std::tie(W0,W1)=BcWeightFactors(v);
      jac+=JacCol(v.mR(),W0*eye(v));
      jac+=JacCol(v.right()->mR(),W1*eye(v));
    }
    else{
      d WR1,WR2; std::tie(WR1,WR2)=BcWeightFactors(v);
      jac+=JacCol(v.mL(),WR1*eye(v));
      jac+=JacCol(v.left()->mL(),WR2*eye(v));
    }
    return jac;
  }

} // Namespace tube

