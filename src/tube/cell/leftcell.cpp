#include "leftcell.h"
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
      d vxi=t[0].vx;
      d vxip1=t[1].vx;

      // Compute weight factors
      d wL0=vxip1/(vxip1-vxi);
      d wL1=-vxi/(vxip1-vxi);
      // VARTRACE(25,wL0);
      // VARTRACE(25,wL1);

      return vd2({wL0,wL1});
    }
  
    vd extrapolateMomentumFlow(const Tube& t){
      vd2 w=weightfactors(t); d wL0=w(0),wL1=w(1);
      return wL0*t[0].mu()()+wL1*t[1].mu()();
    }
    JacRow dExtrapolateMomentumFlow(const Tube& t){
      vd2 w=weightfactors(t); d wL0=w(0),wL1=w(1);
      JacRow jacrow(2);
      jacrow+=JacCol(t[0].mu(),wL0*eye);
      jacrow+=JacCol(t[1].mu(),wL1*eye);
      return jacrow;
    }
  
  } // namespace 
  // Extrapolate momentum flow left side


  
  LeftCell::LeftCell(us i,const Tube& t):
    BcCell(i,t)
  {
    TRACE(15,"LeftCell::LeftCell()");
  }
  void LeftCell::init(const Cell* left,const Cell* right){
    TRACE(10,"LeftCell::init()");
    Cell::init(left,right);
    assert(!left);
    assert(right);

    TL_=var(*gc);
    TL_.setadata(0,gc->T0());
    vars.push_back(&TL_);
    // Remove momentum equation from list. Put equation for Mu in
    // place of momentum eq.
    // WARN("HERE SOME EQS NEED TO BE DELETED");
    assert(eqs.find(EqType::Mom)!=eqs.end());
    // If these elements are already deleted, we do nothing
    try{
      Equation* mom=eqs.at(EqType::Mom);
      delete mom;
      TRACE(25,"MOm eq deleted");
      eqs.erase(EqType::Mom);
    }
    catch(std::out_of_range&){}
    try{
      Equation* mH=eqs.at(EqType::mH_is_m_H);
      delete mH;
      TRACE(25,"mH eq deleted");
      eqs.erase(EqType::mH_is_m_H);
    }
    catch(std::out_of_range&){}
  }
  void LeftCell::show(us detailnr) const{
    cout << "------------- LeftCell ----------\n";
    Cell::show(detailnr);
  }
  vd LeftCell::extrapolateQuant(Physquant p) const {
    TRACE(5,"LeftCell::extrapolateQuant()");
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
  JacRow LeftCell::dExtrapolateQuant(Physquant p) const {
    TRACE(5,"LeftCell::dExtrapolateQuant()");
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
        TRACE(25,"Enthalpy");
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

  // void LeftCell::setResVar(Varnr v,cont vd& res){
  //   TRACE(15,"LeftCell::setResVar()");
    // switch(v){
    // case Varnr::rhoL:
    //   rhoL_.set(res);
    //   break;
    // case Varnr::TL:
    //   TL_.set(res);
    //   break;
    // case Varnr::UL:
    //   UL_.set(res);
    //   break;
    // case Varnr::TsL:
    //   TsL_.set(res);
    //   break;
    // case Varnr::pL:
    //   pL_.set(res);
    //   break;
    // default:
    //   Cell::setResVar(v,res);
    //   break;
    // }
  // }
  // d LeftCell::getValueBc(Varnr v,us freqnr) const{
  //   TRACE(15,"LeftCell::getValueBc()");
    // switch(v) {
    // case Varnr::rho: // Density
    //   return rhoL()(freqnr);
    //   break;
    // case Varnr::U:                 // Volume flown
    //   return UL()(freqnr);
    //   break;
    // case Varnr::p:                   // Pressure
    //   return pL()(freqnr);
    //   break;
    // case Varnr::T:                 // Temp
    //   return TL()(freqnr);
    //   break;
    // case Varnr::Ts:                 // Temp
    //   return TsL()(freqnr);
    //   break;
    // default:
    //   WARN("Unknown variable!");
    //   return 0;
    // }
  // }

  // vd LeftCell::extrapolateMassFlow() const{
  //   TRACE(15,"LeftCell::extrapolateMassFlow()");
  //   const WeightFactors& w=weightFactors();
  //   return w.wL1*right()->continuity().massFlow()+  \
  //    w.wL0*c.massFlow();
  // }
  // JacRow LeftCell::dExtrapolateMassFlow() const{
  //   TRACE(15,"LeftCell::dExtrapolateMassFlow()");
  //   const WeightFactors& w=weightFactors();
  //   JacRow jacrow(-1,4);
  //   jacrow+=JacCol(U(),w.wL0*fDFT*rho().diagt()*iDFT);
  //   jacrow+=JacCol(rho(),w.wL0*fDFT*U().diagt()*iDFT);
  //   jacrow+=JacCol(UR(),w.wL1*fDFT*rhoR().diagt()*iDFT);
  //   jacrow+=JacCol(rhoR(),w.wL1*fDFT*UR().diagt()*iDFT);
  //   return jacrow;
  // }
  // vd LeftCell::extrapolateRhoRT() const{
  //   TRACE(15,"LeftCell::extrapolateRhoRT()");
  //   const WeightFactors& w=weightFactors();
  //   const d& R=gc->gas().Rs();
  //   return fDFT*(R*(w.wL1*rhoR().tdata()%TR().tdata()+  \
  //                              w.wL0*rho().tdata()%T().tdata()));
  // }
  // JacRow LeftCell::dExtrapolateRhoRT() const{
  //   TRACE(15,"LeftCell::dExtrapolateRhoRT()");

  //   const WeightFactors& w=weightFactors();
  //   JacRow jacrow(-1,4);
  //   const d& R=gc->gas().Rs();
  //   jacrow+=JacCol(rho(),fDFT*diagmat(R*w.wL0*T().tdata())*iDFT);
  //   jacrow+=JacCol(T(),fDFT*diagmat(R*w.wL0*rho().tdata())*iDFT);
  //   jacrow+=JacCol(TR(),fDFT*diagmat(R*w.wL1*rhoR().tdata())*iDFT);
  //   jacrow+=JacCol(rhoR(),fDFT*diagmat(R*w.wL1*TR().tdata())*iDFT);
  //   return jacrow;
  // }
  // vd LeftCell::extrapolateMomentumFlow() const{
  //   TRACE(15,"LeftCell::extrapolateMomentumFlow()");
  //   const WeightFactors& w=weightFactors();
  //   return w.wL1*right()->momentum().momentumFlow()+\
  //     w.wL0*m.momentumFlow();
  // }
  // JacRow LeftCell::dExtrapolateMomentumFlow() const{
  //   TRACE(15,"LeftCell::dExtrapolateMomentumFlow()");
  //   const WeightFactors& w=weightFactors();
  //   JacRow jacrow(-1,4);
  //   dmat Utd=U().diagt();
  //   dmat rhotd=rho().diagt();

  //   jacrow+=JacCol(U(),2.0*(w.wL0/w.vSf)*fDFT  \
  //             *rhotd*Utd*iDFT);
  //   jacrow+=JacCol(rho(),(w.wL0/w.vSf)*fDFT\
  //               *(Utd%Utd)*iDFT);

  //   dmat UtdR=UR().diagt();
  //   dmat rhotdR=rhoR().diagt();
  //   jacrow+=JacCol(UR(),2.0*(w.wL1/w.vSfR)*fDFT\
  //              *rhotdR*UtdR*iDFT);
  //   jacrow+=JacCol(rhoR(),(w.wL1/w.vSfR)*fDFT\
  //                *UtdR*UtdR*iDFT);

  //   return jacrow;
  // }

} // namespace tube
