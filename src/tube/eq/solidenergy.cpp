#include "solidenergy.h"
#include "jacobian.h"
#include "cell.h"

#define iDFT (v.gc->iDFT)
#define fDFT (v.gc->fDFT)


namespace tube{
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;


  void SolidTPrescribed::init() {
    TRACE(6,"SolidTPrescribed::init(t)");

  }
  JacRow SolidTPrescribed::jac() const{
    TRACE(6,"SolidTPrescribed::jac()");
    JacRow jac(dofnr,1);
    TRACE(0,"Dofnr jac:"<< dofnr);
    jac+=dTsi();
    return jac;
  }
  vd SolidTPrescribed::error() const {		// Error in momentum equation
    TRACE(6,"SolidTPrescribed::Error()");
    vd error(v.gc->Ns(),fillwith::zeros);
    assert(v.gc!=nullptr);
    error=v.Ts()();
    error(0)-=Tsmirror;
    return error;
  }
  vd SolidTPrescribed::extrapolateHeatFlow() const{
    TRACE(4,"SolidTPrescribed::extrapolateHeatFlow()");
    // const WeightFactors& w=v.weightFactors();
    // vd Qb(v.gc->Ns());
    // if(!v.left()){
    //   vd kappaLt=kappaL();
    //   d SsL;
    //   if(w.SsL>0)
    //     SsL=w.SsL;
    //   else
    //     SsL=1.0;
    //   Qb=(SsL/w.vx)*fDFT*(kappaLt%(v.TsL().tdata()-v.Ts().tdata()));
    // }
    // else if(!v.right()){
    //   vd kappaRt=kappaR();
    //   d SsR;
    //   if(w.SsR>0)
    //     SsR=w.SsR;
    //   else
    //     SsR=1.0;
    //   Qb=(SsR/(w.xR-w.vx))*fDFT*(kappaRt%(v.T().tdata()-v.TR().tdata()));
    // }
    // else{
    //   WARN("That went fatally wrong!");
    //   abort();
    // }
    // return Qb;
  }
  JacRow SolidTPrescribed::dExtrapolateHeatFlow() const{
    TRACE(4,"SolidTPrescribed::dExtrapolateHeatFlow()");
    // JacRow dQb(-1,2);

    // if(!v.left()){
    //   d SsL;
    //   if(w.SsL>0)
    //     SsL=w.SsL;
    //   else
    //     SsL=1.0;
    //   vd kappaLt=kappaL();
    //   dQb+=JacCol(v.Ts(),-(SsL/w.vx)*fDFT*diagmat(kappaLt)*iDFT);
    //   dQb+=JacCol(v.TsL(),(SsL/w.vx)*fDFT*diagmat(kappaLt)*iDFT);
    // }
    // else if(!v.right()){
    //   vd kappaRt=kappaR();
    //   d SsR;
    //   if(w.SsR>0)
    //     SsR=w.SsR;
    //   else
    //     SsR=1.0;
    //   dQb+=JacCol(v.Ts(),(SsR/(w.xR-w.vx))*fDFT*diagmat(kappaRt)*iDFT);
    //   dQb+=JacCol(v.TsR(),-(SsR/(w.xR-w.vx))*fDFT*diagmat(kappaRt)*iDFT);
    // }
    // else{
    //   WARN("That went fatally wrong!");
    //   abort();
    // }
    // return dQb;
  }

  JacCol SolidTPrescribed::dTsi() const {
    TRACE(0,"SolidTPrescribed:dTsi()");
    // Set solid temperature to zero
    return JacCol(v.Ts(),arma::eye(v.gc->Ns(),v.gc->Ns()));
  }
  vd SolidTPrescribed::kappaL() const{
    return ones(v.gc->Ns());
  }
  vd SolidTPrescribed::kappaR() const{
    return ones(v.gc->Ns());
  }

} // namespace tube


