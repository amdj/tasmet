#include "leftcell.h"
#include "weightfactors.h"
#include "jacobian.h"

#define iDFT (gc->iDFT)
#define fDFT (gc->fDFT)

namespace tube{
  using variable::var;
  using tasystem::JacRow;
  using tasystem::JacCol;

  LeftCell::LeftCell(us i,const Tube& t):
    TubeBcCell(i,t)
  {
    TRACE(15,"LeftCell::LeftCell()");
    const tasystem::Globalconf& gc=*(this->gc);
    rhoL_=var(gc);
    UL_=var(gc);
    TL_=var(gc);
    TsL_=var(gc);
  }
  void LeftCell::init(const Cell* left,const Cell* right){
    TRACE(10,"LeftCell::init()");
    Cell::init(left,right);
    assert(!left);
    assert(right);


    vars.push_back(&rhoL_);
    vars.push_back(&UL_);
    vars.push_back(&TL_);
    vars.push_back(&TsL_);

    // Initialize left wall variables
    rhoL_=rho_;
    TL_=T_;
    TsL_=Ts_;
  }
  void LeftCell::show(us detailnr) const{
    cout << "------------- LeftCell ----------\n";
    Cell::show(detailnr);
  }
  void LeftCell::setResVar(varnr v,const vd& res){
    TRACE(15,"LeftCell::setResVar()");
    switch(v){
    case varnr::rhoL:
      rhoL_.set(res);
      break;
    case varnr::TL:
      TL_.set(res);
      break;
    case varnr::UL:
      UL_.set(res);
      break;
    case varnr::TsL:
      TsL_.set(res);
      break;
    case varnr::pL:
      pL_.set(res);
      break;
    default:
      Cell::setResVar(v,res);
      break;
    }
  }
  d LeftCell::getValueBc(varnr v,us freqnr) const{
    TRACE(15,"LeftCell::getValueBc()");
    switch(v) {
    case varnr::rho: // Density
      return rhoL()(freqnr);
      break;
    case varnr::U:                 // Volume flown
      return UL()(freqnr);
      break;
    case varnr::p:                   // Pressure
      return pL()(freqnr);
      break;
    case varnr::T:                 // Temp
      return TL()(freqnr);
      break;
    case varnr::Ts:                 // Temp
      return TsL()(freqnr);
      break;
    default:
      WARN("Unknown variable!");
      return 0;
    }
  }

  vd LeftCell::extrapolateMassFlow() const{
    TRACE(15,"LeftCell::extrapolateMassFlow()");
    const WeightFactors& w=weightFactors();
    return w.wL1*right()->continuity().massFlow()+  \
     w.wL0*c.massFlow();
  }
  JacRow LeftCell::dExtrapolateMassFlow() const{
    TRACE(15,"LeftCell::dExtrapolateMassFlow()");
    const WeightFactors& w=weightFactors();
    JacRow jacrow(-1,4);
    jacrow+=JacCol(U(),w.wL0*fDFT*rho().diagt()*iDFT);
    jacrow+=JacCol(rho(),w.wL0*fDFT*U().diagt()*iDFT);
    jacrow+=JacCol(UR(),w.wL1*fDFT*rhoR().diagt()*iDFT);
    jacrow+=JacCol(rhoR(),w.wL1*fDFT*UR().diagt()*iDFT);
    return jacrow;
  }
  vd LeftCell::extrapolateRhoRT() const{
    TRACE(15,"LeftCell::extrapolateRhoRT()");
    const WeightFactors& w=weightFactors();
    const d& R=gc->gas().Rs();
    return fDFT*(R*(w.wL1*rhoR().tdata()%TR().tdata()+  \
                               w.wL0*rho().tdata()%T().tdata()));
  }
  JacRow LeftCell::dExtrapolateRhoRT() const{
    TRACE(15,"LeftCell::dExtrapolateRhoRT()");

    const WeightFactors& w=weightFactors();
    JacRow jacrow(-1,4);
    const d& R=gc->gas().Rs();
    jacrow+=JacCol(rho(),fDFT*diagmat(R*w.wL0*T().tdata())*iDFT);
    jacrow+=JacCol(T(),fDFT*diagmat(R*w.wL0*rho().tdata())*iDFT);
    jacrow+=JacCol(TR(),fDFT*diagmat(R*w.wL1*rhoR().tdata())*iDFT);
    jacrow+=JacCol(rhoR(),fDFT*diagmat(R*w.wL1*TR().tdata())*iDFT);
    return jacrow;
  }
  vd LeftCell::extrapolateMomentumFlow() const{
    TRACE(15,"LeftCell::extrapolateMomentumFlow()");
    const WeightFactors& w=weightFactors();
    return w.wL1*right()->momentum().momentumFlow()+\
      w.wL0*m.momentumFlow();
  }
  JacRow LeftCell::dExtrapolateMomentumFlow() const{
    TRACE(15,"LeftCell::dExtrapolateMomentumFlow()");
    const WeightFactors& w=weightFactors();
    JacRow jacrow(-1,4);
    dmat Utd=U().diagt();
    dmat rhotd=rho().diagt();

    jacrow+=JacCol(U(),2.0*(w.wL0/w.vSf)*fDFT  \
              *rhotd*Utd*iDFT);
    jacrow+=JacCol(rho(),(w.wL0/w.vSf)*fDFT\
                *(Utd%Utd)*iDFT);

    dmat UtdR=UR().diagt();
    dmat rhotdR=rhoR().diagt();
    jacrow+=JacCol(UR(),2.0*(w.wL1/w.vSfR)*fDFT\
               *rhotdR*UtdR*iDFT);
    jacrow+=JacCol(rhoR(),(w.wL1/w.vSfR)*fDFT\
                 *UtdR*UtdR*iDFT);

    return jacrow;
  }

} // namespace tube
