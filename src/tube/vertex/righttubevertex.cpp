#include "righttubevertex.h"
#include "weightfactors.h"
#include "jacobian.h"

#define iDFT (gc->iDFT)
#define fDFT (gc->fDFT)


namespace tube{
  using variable::var;
  using tasystem::JacRow;
  using tasystem::JacCol;

  RightTubeVertex::RightTubeVertex(us i,const Tube& t):
    TubeBcVertex(i,t),
    sR(*this)
  {
    // Initialize right wall variables
    rhoR_=var(gc);
    UR_=var(gc);
    TR_=var(gc);
    pR_=var(gc);
    TsR_=var(gc);
  }
  void RightTubeVertex::init(const TubeVertex* left,const TubeVertex* right){
    TRACE(15,"RightTubeVertex::init()");
    assert(!right);             // Otherwise, this is not the
                                // rightmost!
    assert(left);
    TubeVertex::init(left,right);
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
  void RightTubeVertex::show(us detailnr) const{
    cout << "------------- RightTubeVertex ---------\n";
    TubeVertex::show(detailnr);
  }
  void RightTubeVertex::setResVar(varnr v,const vd& res){
    TRACE(15,"RightTubeVertex::setResVar()");
    switch(v){
    case varnr::rhoR:
      rhoR_.set(res);
      break;
    case varnr::TR:
      TR_.set(res);
      break;
    case varnr::UR:
      UR_.set(res);
      break;
    case varnr::TsR:
      TsR_.set(res);
      break;
    case varnr::pR:
      pR_.set(res);
      break;
    default:
      TubeVertex::setResVar(v,res);
      break;
    }
  }
  vd RightTubeVertex::extrapolateMassFlow() const{
    TRACE(15,"RightTubeVertex::extrapolateMassFlow()");
    const WeightFactors& w=weightFactors();
    return w.wRNm2*left()->continuity().massFlow()+\
      w.wRNm1*c.massFlow();
  }
  JacRow RightTubeVertex::dExtrapolateMassFlow() const{
    TRACE(15,"RightTubeVertex::dExtrapolateMassFlow()");
    const WeightFactors& w=weightFactors();
    JacRow jacrow(-1,4);
    jacrow+=JacCol(U(),w.wRNm1*fDFT*rho().diagt()*iDFT);
    jacrow+=JacCol(rho(),w.wRNm1*fDFT*U().diagt()*iDFT);
    jacrow+=JacCol(UL(),w.wRNm2*fDFT*rhoL().diagt()*iDFT);
    jacrow+=JacCol(rhoL(),w.wRNm2*fDFT*UL().diagt()*iDFT);
    return jacrow;
  }
  vd RightTubeVertex::extrapolateRhoRT() const{
    const WeightFactors& w=weightFactors();
    const d& R=gc->gas.Rs();
    return fDFT*(R*(w.wRNm2*rhoL().tdata()%TL().tdata()+   \
                    w.wRNm1*rho().tdata()%T().tdata()));
  }
  JacRow RightTubeVertex::dExtrapolateRhoRT() const{
    TRACE(15,"RightTubeVertex::dExtrapolateRhoRT()");
    const WeightFactors& w=weightFactors();
    JacRow jacrow(-1,4);
    const d& R=gc->gas.Rs();    
    jacrow+=JacCol(rho(),R*w.wRNm1*fDFT*diagmat(T().tdata())*iDFT);
    jacrow+=JacCol(T(),R*w.wRNm1*fDFT*diagmat(rho().tdata())*iDFT);
    jacrow+=JacCol(rhoL(),R*w.wRNm2*fDFT*diagmat(TL().tdata())*iDFT);
    jacrow+=JacCol(TL(),R*w.wRNm2*fDFT*diagmat(rhoL().tdata())*iDFT);
    return jacrow;
  }
  vd RightTubeVertex::extrapolateMomentumFlow() const{
    TRACE(15,"RightTubeVertex::extrapolateMomentumFlow()");
    const WeightFactors& w=weightFactors();
    return w.wRNm2*left()->momentum().momentumFlow()+\
      w.wRNm1*m.momentumFlow();
  }
  JacRow RightTubeVertex::dExtrapolateMomentumFlow() const{
    TRACE(15,"RightTubeVertex::dExtrapolateMomentumFlow()");
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
  d RightTubeVertex::getResBc(varnr v,us freqnr) const{
    TRACE(15,"RightTubeVertex::getResBc()");
    switch(v) {
    case varnr::rho: // Density
      return rhoR()(freqnr);
      break;
    case varnr::U:                 // Volume flown
      return UR()(freqnr);
      break;
    case varnr::p:                   // Pressure
      return pR()(freqnr);
      break;
    case varnr::T:                 // Temp
      return TR()(freqnr);
      break;
    case varnr::Ts:                 // Temp
      return TsR()(freqnr);
      break;
    default:
      WARN("Unknown variable!");
      return 0;
    }
  }

} // namespace tube

