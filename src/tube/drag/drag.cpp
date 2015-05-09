#include "drag.h"
#include "cell.h"



namespace tube{
  namespace drag {

    using std::make_tuple;
    using std::tuple;

    tuple<d,d> wf(const Cell& v) {
      TRACE(5,"wf(Cell)");
      assert(v.left());
      const d& xim1=v.left()->vx;
      const d& xi=v.vx;
      const d& xL=v.xL;
      const d wim1=(xi-xL)/(xi-xim1);
      return make_tuple(wim1,1-wim1);
    }

    vd DragResistance::drag(const Cell& v) const {return vd(v.gc->Ns(),fillwith::zeros);}
    dmat DragResistance::dm(const Cell& v) const {return zeros<dmat>(v.gc->Ns(),v.gc->Ns());}
    // dmat DragResistance::drhoi(const Cell& v) const {return v.zero;}
    // dmat DragResistance::dpi(const Cell& v) const {return v.zero;}
    vd DragResistance::shearWaveNumber(const Cell& v) const {
      TRACE(15,"DragResistance::shearWaveNumber(const Cell& v)");
      const us Nf=v.gc->Nf();
      d T0,p0;
      if(v.left()) {
        d wim1,wi; std::tie(wim1,wi)=wf(v);
        T0=wi*v.T()(0)+wim1*v.TL()(0);	// Time-averaged temperature at the cell wall
        p0=v.gc->p0()+wi*v.p()(0)+wim1*v.p()(0);	// Time-averaged temperature at the cell wall
      }
      else {
        T0=v.T()(0);	// Time-averaged temperature at the cell wall
        p0=v.gc->p0()+v.p()(0);	// Time-averaged temperature at the cell wall
      }
      const d mu0=v.gc->gas().mu(T0);
      const d rho0=v.gc->gas().rho(T0,p0);

      const d omg=v.gc->getomg();
      const vd omgvec=linspace(0,Nf*omg,Nf+1);
      TRACE(20,"SFSQ")
      const d& rh=v.rhL;

      return rh*sqrt((rho0/mu0)*omgvec);

    }
  } // namespace drag
} // namespace tube





















