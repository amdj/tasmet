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
  vd RightTubeVertex::extrapolateDensity() const{
    const WeightFactors& w=weightFactors();
    return w.wRNm2*rhoL()()+                    \
      w.wRNm1*rho()();
  }
  // vd RightTubeVertex::extrapolateVolumeFlow() const{
  //   const WeightFactors& w=weightFactors();
  //   return w.wRNm2*UL()()+                    \
  //     w.wRNm1*U()();
  // }
  JacRow RightTubeVertex::dExtrapolateDensity() const{
    TRACE(15,"RightTubeVertex::dExtrapolateDensity()");
    const WeightFactors& w=weightFactors();
    JacRow jacrow(-1,2);
    jacrow+=JacCol(rho(),w.wRNm1*eye<dmat>(gc->Ns(),gc->Ns()));
    jacrow+=JacCol(rhoL(),w.wRNm2*eye<dmat>(gc->Ns(),gc->Ns()));
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

} // namespace tube

