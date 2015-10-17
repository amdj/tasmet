#include "cell.h"
#include "jacrow.h"
#include "var.h"
#include "continuity.h"
#include "weightfactors.h"

#define iDFT (v.gc->iDFT())
#define fDFT (v.gc->fDFT())
#define DDTfd (v.gc->DDTfd())
#define Ns (v.gc->Ns())
namespace duct{

  using tasystem::JacCol;
  using tasystem::JacRow;
  using tasystem::var;
  
  void Continuity::show() const {
    cout << "-------------- Continuity equation\n";
    cout << "Wddt : " << Wddt << "\n";
  }
  Continuity::~Continuity(){TRACE(0,"~Continuity()");}
  void Continuity::init(){
    TRACE(8,"Continuity::init(duct)");
    Wddt=v.vVf;
  } // init
  JacRow Continuity::jac() const{
    TRACE(15,"Continuity::jac()");
    JacRow jac(dofnr,3);
    // TRACE(0,"Continuity, dofnr jac:"<< dofnr);
    jac+=JacCol(v.rho(),Wddt*DDTfd);    
    jac+=JacCol(v.mr(),eye());
    jac+=JacCol(v.ml(),-eye());
    return jac;
  }
  vd Continuity::error() const {	
    TRACE(15,"Continuity::Error()");
    vd error(Ns,fillwith::zeros);
    error+=Wddt*DDTfd*v.rho()();
    error+=v.mr()();
    error-=v.ml()();

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
      d W0,W1; std::tie(W0,W1)=BcWeightFactorsW(v);
      return W0*v.mr()()+W1*v.right()->mr()();
    }
    else{
      d WR1,WR2; std::tie(WR1,WR2)=BcWeightFactorsW(v);
      return WR1*v.ml()()+WR2*v.left()->ml()();
    }
  }
  JacRow Continuity::dExtrapolateMassFlow(const Cell& v){
    TRACE(15,"Continuity::dExtrapolateMassFlow(const Cell& v)");
    assert((!v.left() && v.right()) || (!v.right() && v.left()));
    JacRow jac(2);
    if(!v.left()){
      d W0,W1; std::tie(W0,W1)=BcWeightFactorsW(v);
      jac+=JacCol(v.mr(),W0*eye(v));
      jac+=JacCol(v.right()->mr(),W1*eye(v));
    }
    else{
      d WR1,WR2; std::tie(WR1,WR2)=BcWeightFactorsW(v);
      jac+=JacCol(v.ml(),WR1*eye(v));
      jac+=JacCol(v.left()->ml(),WR2*eye(v));
    }
    return jac;
  }

} // Namespace duct

