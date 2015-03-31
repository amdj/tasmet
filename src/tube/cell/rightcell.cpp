#include "rightcell.h"
#include "weightfactors.h"
#include "jacobian.h"

#define iDFT (gc->iDFT)
#define fDFT (gc->fDFT)


namespace tube{
  using variable::var;
  using tasystem::JacRow;
  using tasystem::JacCol;

  RightCell::RightCell(us i,const Tube& t):
    BcCell(i,t),
    sR(*this)
  {
    const tasystem::Globalconf& gc=*(this->gc);
    // Initialize right wall variables
    rhoR_=var(gc);
    UR_=var(gc);
    TR_=var(gc);
    pR_=var(gc);
    TsR_=var(gc);
  }
  void RightCell::init(const Cell* left,const Cell* right){
    TRACE(15,"RightCell::init()");
    assert(!right);             // Otherwise, this is not the
                                // rightmost!
    assert(left);
    Cell::init(left,right);
    sR.init();

    eqs.push_back(&sR);

    vars.push_back(&rhoR_);
    vars.push_back(&UR_);
    vars.push_back(&TR_);
    vars.push_back(&TsR_);
    vars.push_back(&pR_);    
    pR_=pL_;
    rhoR_=rho_;
    TR_=T_;
    TsR_=Ts_;
  }
  void RightCell::show(us detailnr) const{
    cout << "------------- RightCell ---------\n";
    Cell::show(detailnr);
  }
  void RightCell::setResVar(Varnr v,const vd& res){
    TRACE(15,"RightCell::setResVar()");
    switch(v){
    case Varnr::rhoR:
      rhoR_.set(res);
      break;
    case Varnr::TR:
      TR_.set(res);
      break;
    case Varnr::UR:
      UR_.set(res);
      break;
    case Varnr::TsR:
      TsR_.set(res);
      break;
    case Varnr::pR:
      pR_.set(res);
      break;
    default:
      Cell::setResVar(v,res);
      break;
    }
  }
  vd RightCell::extrapolateMassFlow() const{
    TRACE(15,"RightCell::extrapolateMassFlow()");
    const WeightFactors& w=weightFactors();
    return w.wRNm2*left()->continuity().massFlow()+\
      w.wRNm1*c.massFlow();
  }
  JacRow RightCell::dExtrapolateMassFlow() const{
    TRACE(15,"RightCell::dExtrapolateMassFlow()");
    const WeightFactors& w=weightFactors();
    JacRow jacrow(-1,4);
    jacrow+=JacCol(U(),w.wRNm1*fDFT*rho().diagt()*iDFT);
    jacrow+=JacCol(rho(),w.wRNm1*fDFT*U().diagt()*iDFT);
    jacrow+=JacCol(UL(),w.wRNm2*fDFT*rhoL().diagt()*iDFT);
    jacrow+=JacCol(rhoL(),w.wRNm2*fDFT*UL().diagt()*iDFT);
    return jacrow;
  }
  vd RightCell::extrapolateRhoRT() const{
    const WeightFactors& w=weightFactors();
    const d& R=gc->gas().Rs();
    return fDFT*(R*(w.wRNm2*rhoL().tdata()%TL().tdata()+   \
                    w.wRNm1*rho().tdata()%T().tdata()));
  }
  JacRow RightCell::dExtrapolateRhoRT() const{
    TRACE(15,"RightCell::dExtrapolateRhoRT()");
    const WeightFactors& w=weightFactors();
    JacRow jacrow(-1,4);
    const d& R=gc->gas().Rs();    
    jacrow+=JacCol(rho(),R*w.wRNm1*fDFT*diagmat(T().tdata())*iDFT);
    jacrow+=JacCol(T(),R*w.wRNm1*fDFT*diagmat(rho().tdata())*iDFT);
    jacrow+=JacCol(rhoL(),R*w.wRNm2*fDFT*diagmat(TL().tdata())*iDFT);
    jacrow+=JacCol(TL(),R*w.wRNm2*fDFT*diagmat(rhoL().tdata())*iDFT);
    return jacrow;
  }
  vd RightCell::extrapolateMomentumFlow() const{
    TRACE(15,"RightCell::extrapolateMomentumFlow()");
    const WeightFactors& w=weightFactors();
    return w.wRNm2*left()->momentum().momentumFlow()+\
      w.wRNm1*m.momentumFlow();
  }
  JacRow RightCell::dExtrapolateMomentumFlow() const{
    TRACE(15,"RightCell::dExtrapolateMomentumFlow()");
    const WeightFactors& w=weightFactors();
    JacRow jacrow(-1,4);
    dmat Utd=U().diagt();
    dmat rhotd=rho().diagt();

    jacrow+=JacCol(U(),2.0*w.wRNm1/w.vSf*fDFT\
                   *rhotd*Utd*iDFT);
    jacrow+=JacCol(rho(),(w.wRNm1/w.vSf)*fDFT\
                   *(Utd%Utd)*iDFT);

    dmat UtdL=UL().diagt();
    dmat rhotdL=rhoL().diagt();
    jacrow+=JacCol(UL(),2.0*(w.wRNm2/w.vSfL)*fDFT\
                   *rhotdL*UtdL*iDFT);
    jacrow+=JacCol(rhoL(),(w.wRNm2/w.vSfL)*fDFT\
                   *UtdL*UtdL*iDFT);
    return jacrow;
  }
  d RightCell::getValueBc(Varnr v,us freqnr) const{
    TRACE(15,"RightCell::getValueBc()");
    switch(v) {
    case Varnr::rho: // Density
      return rhoR()(freqnr);
      break;
    case Varnr::U:                 // Volume flown
      return UR()(freqnr);
      break;
    case Varnr::p:                   // Pressure
      return pR()(freqnr);
      break;
    case Varnr::T:                 // Temp
      return TR()(freqnr);
      break;
    case Varnr::Ts:                 // Temp
      return TsR()(freqnr);
      break;
    default:
      WARN("Unknown variable!");
      return 0;
    }
  }

} // namespace tube

