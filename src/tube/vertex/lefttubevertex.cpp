#include "lefttubevertex.h"
#include "weightfactors.h"
#include "jacobian.h"

#define iDFT (gc->iDFT)
#define fDFT (gc->fDFT)

namespace tube{
  using variable::var;
  using tasystem::JacRow;
  using tasystem::JacCol;

  LeftTubeVertex::LeftTubeVertex(us i,const Tube& t):
    TubeBcVertex(i,t)
  {
    TRACE(15,"LeftTubeVertex::LeftTubeVertex()");
    rhoL_=var(gc);
    UL_=var(gc);
    TL_=var(gc);
    TsL_=var(gc);
  }
  void LeftTubeVertex::init(const TubeVertex* left,const TubeVertex* right){
    TRACE(10,"LeftTubeVertex::init()");
    TubeVertex::init(left,right);
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
  void LeftTubeVertex::show(us detailnr) const{
    cout << "------------- LeftTubeVertex ----------\n";
    TubeVertex::show(detailnr);
  }
  void LeftTubeVertex::setResVar(varnr v,const vd& res){
    TRACE(15,"LeftTubeVertex::setResVar()");
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
      TubeVertex::setResVar(v,res);
      break;
    }
  }
  vd LeftTubeVertex::extrapolateMassFlow() const{
    TRACE(15,"LeftTubeVertex::extrapolateMassFlow()");
    const WeightFactors& w=weightFactors();
    return w.wL1*right()->continuity().massFlow()+  \
     w.wL0*c.massFlow();
  }
  JacRow LeftTubeVertex::dExtrapolateMassFlow() const{
    TRACE(15,"LeftTubeVertex::dExtrapolateMassFlow()");
    const WeightFactors& w=weightFactors();
    JacRow jacrow(-1,4);
    jacrow+=JacCol(U(),w.wL0*fDFT*rho().diagt()*iDFT);
    jacrow+=JacCol(rho(),w.wL0*fDFT*U().diagt()*iDFT);
    jacrow+=JacCol(UR(),w.wL1*fDFT*rhoR().diagt()*iDFT);
    jacrow+=JacCol(rhoR(),w.wL1*fDFT*UR().diagt()*iDFT);
    return jacrow;
  }
  vd LeftTubeVertex::extrapolateDensity() const{
    TRACE(15,"LeftTubeVertex::extrapolateDensity()");
    const WeightFactors& w=weightFactors();
    return w.wL1*rhoR()()+\
      w.wL0*rho()();
  }
  JacRow LeftTubeVertex::dExtrapolateDensity() const{
    TRACE(15,"LeftTubeVertex::dExtrapolateDensity()");

    const WeightFactors& w=weightFactors();
    JacRow jacrow(-1,2);
    jacrow+=JacCol(rho(),w.wL0*eye<dmat>(gc->Ns(),gc->Ns()));
    jacrow+=JacCol(rhoR(),w.wL1*eye<dmat>(gc->Ns(),gc->Ns()));
    return jacrow;
  }
  vd LeftTubeVertex::extrapolateMomentumFlow() const{
    TRACE(15,"LeftTubeVertex::extrapolateMomentumFlow()");
    const WeightFactors& w=weightFactors();
    return w.wL1*right()->momentum().momentumFlow()+\
      w.wL0*m.momentumFlow();
  }
  JacRow LeftTubeVertex::dExtrapolateMomentumFlow() const{
    TRACE(15,"LeftTubeVertex::dExtrapolateMomentumFlow()");
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
