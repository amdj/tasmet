#include "bccell.h"
#include "continuity.h"
#include "momentum.h"
#include "energy.h"
#include "state.h"

#include "jacrow.h"
#include "exception.h"
#define iDFT (v.gc->iDFT)
#define fDFT (v.gc->fDFT)
#define DDTfd (v.gc->DDTfd)
#define Ns (v.gc->Ns())

namespace tube{

  using tasystem::JacRow;
  using tasystem::JacCol;
  using tasystem::var;

  BcCell::BcCell(us i,const Tube& t):
    Cell(i,t)
  {}
  void BcCell::init(const Cell* left,const Cell* right){
    TRACE(10,"BcCell::init(const Cell* left,const Cell* right)");
    rhobc_=var(*gc,gc->rho0());
    Tbc_=var(*gc,gc->T0());
    Tsbc_=var(*gc,gc->T0());
    pbc_=var(*gc);
    ubc_=var(*gc);
    mubc_=var(*gc);
    mHbc_=var(*gc);
    Cell::init(left,right);
  }
  vd BcCell::extrapolateQuant(Varnr v) const {
    TRACE(5,"BcCell::extrapolateQuant()");
    switch(v){
    case Varnr::Q:
      return Energy::extrapolateHeatFlow(*this);
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

}                // namespace tube

