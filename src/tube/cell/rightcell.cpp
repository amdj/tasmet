#include "rightcell.h"
#include "jacrow.h"
#include "tube.h"
#include "energy.h"

#define iDFT (gc->iDFT)
#define fDFT (gc->fDFT)
#define eye eye(t.Gc().Ns(),t.Gc().Ns())

namespace tube{
  using variable::var;
  using tasystem::JacRow;
  using tasystem::JacCol;


  // These functions should stay internal to this unit
  namespace {
    
    vd2 weightfactors(const Tube& t){
      us nCells=t.getNCells();
      d xR=t[nCells-1].xR;
      d vxm1=t[nCells-1].vx;
      d vxm2=t[nCells-2].vx;

      // Compute weight factors
      d wRNm1=(vxm2-xR)/(vxm2-vxm1);
      d wRNm2=(xR-vxm1)/(vxm2-vxm1);

      VARTRACE(25,wRNm1);
      VARTRACE(25,wRNm2);
      return vd2({wRNm1,wRNm2});
    }
   // Extrapolate momentum flow left side  
    vd extrapolateMomentumFlow(const Tube& t){
      us nCells=t.getNCells();
      vd2 w=weightfactors(t); d wRNm1=w(0),wRNm2=w(1);
      return wRNm1*t[nCells-1].mu()()+wRNm2*t[nCells-2].mu()();
    }
    JacRow dExtrapolateMomentumFlow(const Tube& t){
      us nCells=t.getNCells();
      vd2 w=weightfactors(t); d wRNm1=w(0),wRNm2=w(1);

      JacRow jacrow(2);
      jacrow+=JacCol(t[nCells-1].mu(),wRNm1*eye);
      jacrow+=JacCol(t[nCells-2].mu(),wRNm2*eye);
      return jacrow;
    }
    
  } // namespace Namespace

 RightCell::RightCell(us i,const Tube& t):
    BcCell(i,t)
  {
    const tasystem::Globalconf& gc=*(this->gc);
    // Initialize right wall variables
  }
  void RightCell::init(const Cell* left,const Cell* right){
    TRACE(15,"RightCell::init()");
    assert(!right);             // Otherwise, this is not the
                                // rightmost!
    assert(left);
    Cell::init(left,right);

    TR_=var(*gc);
    TR_.setadata(0,gc->T0());
    mR_=var(*gc);
    mHR_=var(*gc);
    vars.push_back(&mR_);
    vars.push_back(&mHR_);
    vars.push_back(&TR_);
    
  }
  void RightCell::show(us detailnr) const{
    cout << "------------- RightCell ---------\n";
    Cell::show(detailnr);
  }

  vd RightCell::extrapolateQuant(Physquant p) const {
    TRACE(5,"RightCell::extrapolateQuant()");
    switch(p){
    case Physquant::momentumFlow:
      return extrapolateMomentumFlow(getTube());      
      break;
    case Physquant::heatFlow:
      {
        Energy e(*this);
        e.init();
        return e.extrapolateHeatFlow();
        break;
      }
    case Physquant::enthalpyFlow:
      {
        Energy e(*this);
        e.init();
        return e.extrapolateEnthalpyFlow();
        break;
      }
    default:
      WARN("This is not yet implemented!");
      assert(false);
    }

  }
  JacRow RightCell::dExtrapolateQuant(Physquant p) const {
    TRACE(5,"RightCell::dExtrapolateQuant()");
    switch(p){
    case Physquant::momentumFlow:
      return dExtrapolateMomentumFlow(getTube());      
      break;
    case Physquant::heatFlow:
      {
        Energy e(*this);
        e.init();
        return e.dExtrapolateHeatFlow();
        break;
      }
    case Physquant::enthalpyFlow:
      {
        Energy e(*this);
        e.init();
        return e.dExtrapolateEnthalpyFlow();
        break;
      }
    default:
      WARN("This is not yet implemented!");
      assert(false);

    }
  }

  // void RightCell::setResVar(Varnr v,const vd& res){
  //   TRACE(15,"RightCell::setResVar()");
    // switch(v){
    // case Varnr::rhoR:
    //   rhoR_.set(res);
    //   break;
    // case Varnr::TR:
    //   TR_.set(res);
    //   break;
    // case Varnr::UR:
    //   UR_.set(res);
    //   break;
    // case Varnr::TsR:
    //   TsR_.set(res);
    //   break;
    // case Varnr::pR:
    //   pR_.set(res);
    //   break;
    // default:
    //   Cell::setResVar(v,res);
    //   break;
    // }
  // }
  // vd RightCell::extrapolateMassFlow() const{
  //   TRACE(15,"RightCell::extrapolateMassFlow()");
  //   const WeightFactors& w=weightFactors();
  //   return w.wRNm2*left()->continuity().massFlow()+\
  //     w.wRNm1*c.massFlow();
  // }
  // JacRow RightCell::dExtrapolateMassFlow() const{
  //   TRACE(15,"RightCell::dExtrapolateMassFlow()");
  //   const WeightFactors& w=weightFactors();
  //   JacRow jacrow(-1,4);
  //   jacrow+=JacCol(U(),w.wRNm1*fDFT*rho().diagt()*iDFT);
  //   jacrow+=JacCol(rho(),w.wRNm1*fDFT*U().diagt()*iDFT);
  //   jacrow+=JacCol(UL(),w.wRNm2*fDFT*rhoL().diagt()*iDFT);
  //   jacrow+=JacCol(rhoL(),w.wRNm2*fDFT*UL().diagt()*iDFT);
  //   return jacrow;
  // }
  // vd RightCell::extrapolateRhoRT() const{
  //   const WeightFactors& w=weightFactors();
  //   const d& R=gc->gas().Rs();
  //   return fDFT*(R*(w.wRNm2*rhoL().tdata()%TL().tdata()+   \
  //                   w.wRNm1*rho().tdata()%T().tdata()));
  // }
  // JacRow RightCell::dExtrapolateRhoRT() const{
  //   TRACE(15,"RightCell::dExtrapolateRhoRT()");
  //   const WeightFactors& w=weightFactors();
  //   JacRow jacrow(-1,4);
  //   const d& R=gc->gas().Rs();    
  //   jacrow+=JacCol(rho(),R*w.wRNm1*fDFT*diagmat(T().tdata())*iDFT);
  //   jacrow+=JacCol(T(),R*w.wRNm1*fDFT*diagmat(rho().tdata())*iDFT);
  //   jacrow+=JacCol(rhoL(),R*w.wRNm2*fDFT*diagmat(TL().tdata())*iDFT);
  //   jacrow+=JacCol(TL(),R*w.wRNm2*fDFT*diagmat(rhoL().tdata())*iDFT);
  //   return jacrow;
  // }
  // vd RightCell::extrapolateMomentumFlow() const{
  //   TRACE(15,"RightCell::extrapolateMomentumFlow()");
  //   const WeightFactors& w=weightFactors();
  //   return w.wRNm2*left()->momentum().momentumFlow()+\
  //     w.wRNm1*m.momentumFlow();
  // }
  // JacRow RightCell::dExtrapolateMomentumFlow() const{
  //   TRACE(15,"RightCell::dExtrapolateMomentumFlow()");
  //   const WeightFactors& w=weightFactors();
  //   JacRow jacrow(-1,4);
  //   dmat Utd=U().diagt();
  //   dmat rhotd=rho().diagt();

  //   jacrow+=JacCol(U(),2.0*w.wRNm1/w.vSf*fDFT\
  //                  *rhotd*Utd*iDFT);
  //   jacrow+=JacCol(rho(),(w.wRNm1/w.vSf)*fDFT\
  //                  *(Utd%Utd)*iDFT);

  //   dmat UtdL=UL().diagt();
  //   dmat rhotdL=rhoL().diagt();
  //   jacrow+=JacCol(UL(),2.0*(w.wRNm2/w.vSfL)*fDFT\
  //                  *rhotdL*UtdL*iDFT);
  //   jacrow+=JacCol(rhoL(),(w.wRNm2/w.vSfL)*fDFT\
  //                  *UtdL*UtdL*iDFT);
  //   return jacrow;
  // }
  // d RightCell::getValueBc(Varnr v,us freqnr) const{
  //   TRACE(15,"RightCell::getValueBc()");
  //   switch(v) {
  //   case Varnr::rho: // Density
  //     return rhoR()(freqnr);
  //     break;
  //   case Varnr::U:                 // Volume flown
  //     return UR()(freqnr);
  //     break;
  //   case Varnr::p:                   // Pressure
  //     return pR()(freqnr);
  //     break;
  //   case Varnr::T:                 // Temp
  //     return TR()(freqnr);
  //     break;
  //   case Varnr::Ts:                 // Temp
  //     return TsR()(freqnr);
  //     break;
  //   default:
  //     WARN("Unknown variable!");
  //     return 0;
  //   }
  // }

} // namespace tube


