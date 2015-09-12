#include "bccell.h"
#include "continuity.h"
#include "momentum.h"
#include "energy.h"
#include "solidenergy.h"
#include "state.h"

#include "jacrow.h"
#include "exception.h"
#define iDFT (v.gc->iDFT)
#define fDFT (v.gc->fDFT)
#define DDTfd (v.gc->DDTfd)
#define Ns (gc->Ns())

namespace duct{

  using tasystem::JacRow;
  using tasystem::JacCol;
  using tasystem::var;

  BcCell::BcCell(us i,const Duct& t):
    Cell(i,t)
  {
    rhobc_=var(*gc,gc->rho0());
    Tbc_=var(*gc,gc->T0());
    Tsbc_=var(*gc,gc->T0());
    pbc_=var(*gc);
    ubc_=var(*gc);
    mubc_=var(*gc);
    mHbc_=var(*gc);
}
  void BcCell::init(const Cell* left,const Cell* right){
    TRACE(10,"BcCell::init(const Cell* left,const Cell* right)");
    Cell::init(left,right);    
  }
  void BcCell::updateNf(){
    TRACE(15,"BcCell::updateNf()");
    Cell::updateNf();
    rhobc_.updateNf();
    Tbc_.updateNf();
    Tsbc_.updateNf();
    pbc_.updateNf();
    ubc_.updateNf();
    mubc_.updateNf();
    mHbc_.updateNf();
  }
  vd BcCell::extrapolateQuant(Varnr v) const {
    TRACE(5,"BcCell::extrapolateQuant()");
    switch(v){
    case Varnr::Q:
      return Energy::extrapolateHeatFlow(*this);
      break;
    case Varnr::Qs:
      if(eqs.find(EqType::Sol)!=eqs.end())
	return static_cast<const SolidEnergy*>(eqs.at(EqType::Sol))->extrapolateHeatFlow();
      else
	return vd(Ns,fillwith::zeros);
      break;
    case Varnr::mu:
      return ExtrapolateMomentumFlow::extrapolateMomentumFlow(*this);
      break;
    case Varnr::mH:
      return ExtrapolateEnthalpyFlow::extrapolateEnthalpyFlow(*this);
      break;
    default:
      WARN("This is not implemented!");
      throw MyError("Unextrapolated variable asked");
    }
  }
  JacRow BcCell::dExtrapolateQuant(Varnr p) const {
    TRACE(5,"LeftCell::dExtrapolateQuant()");
    switch(p){
    case Varnr::Q:
      return Energy::dExtrapolateHeatFlow(*this);
      break;
    case Varnr::Qs:
      if(eqs.find(EqType::Sol)!=eqs.end())
	return static_cast<const SolidEnergy*>(eqs.at(EqType::Sol))->dExtrapolateHeatFlow();
      else
	return JacRow(-1,0);
      break;
    case Varnr::mu:
      return ExtrapolateMomentumFlow::dExtrapolateMomentumFlow(*this);
      break;
    case Varnr::mH:
      return ExtrapolateEnthalpyFlow::dExtrapolateEnthalpyFlow(*this);
      break;
    default:
      WARN("This is not implemented!");
      throw MyError("Unextrapolated variable asked");
    }
  }

}                // namespace duct

