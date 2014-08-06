#include "tube.h"
#include "tubevertex.h"
#include "laminardrag.h"

namespace tube{
  
  LaminarDragResistance::LaminarDragResistance(const Tube& t):zfd(t){
    TRACE(10,"LaminarDragResistanc::LaminarDragResistance()");
    rf=rottfuncs::RottFuncs(t.geom.shape); // Reinitialize thermoviscous functions with right shape
  }

  vc LaminarDragResistance::ComplexResistancecoef(const TubeVertex& v) const {
    TRACE(0,"LaminarDragResistance::ComplexResistancecoef()");
    const us& Nf=v.gc->Nf;
    const us& i=v.i;
    d T0=v.T(0);	// Time-averaged temperature
    d mu0=v.gc->gas.mu(T0);
    d p0=v.p(0)+v.gc->p0;
    d rho0=v.gc->gas.rho(T0,p0);

    const d& rh=v.lg.vrh;
    TRACE(-1,"rh: "<<rh);
    d omg,deltanu; 
    vc rh_over_deltanu(Nf);
    vc omgvec(Nf);
    for(us j=0;j<v.gc->Nf;j++)
      {
	TRACE(-1,"j:"<<j);
	omg=v.gc->omg*(j+1);
	omgvec(j)=omg;
	deltanu=sqrt(2*mu0/(rho0*omg));
	TRACE(-1,"deltanu: " << deltanu);
	rh_over_deltanu(j)=rh/deltanu;
      }
    vc fnu=rf.fx(rh_over_deltanu); // Viscous rott function
    TRACE(-1,"fnu:" << fnu);
    TRACE(-1,"omgvec:"<<omgvec);
    vc Resistancecoef=I*rho0*omgvec%(fnu/(1.0-fnu));
    return Resistancecoef;
  }
  vd LaminarDragResistance::drag(const TubeVertex& v) const {
    TRACE(10,"LaminarDragResistance::operator()");
    const us& i=v.i;
    const d& rh=v.lg.vrh;
    const us& Nf=v.gc->Nf;
    dmat dDdU=dUi(v);
    vd drag=dDdU*v.U();
    // VERY IMPORTANT: NOM
    return drag; 		// No momentum scale here, since this is already done in dUi!!!!
  }
  dmat LaminarDragResistance::dUi(const TubeVertex& v) const { // Derivative of drag resistance to velocity
    TRACE(10,"LaminarDragResistance::dUi()");
    dmat dUi(v.gc->Ns,v.gc->Ns,fillwith::zeros);
    const us& i=v.i;
    const d& rh=v.lg.vrh;
    const us& Nf=v.gc->Nf;
    d T0=v.T(0);	// Time-averaged temperature
    d mu0=v.gc->gas.mu(T0);

    // The complex resistance coefficient vector has size Nf
    vc CResistance=ComplexResistancecoef(v);
    for(us j=1;j<Nf+1;j++){
      dUi(2*j-1,2*j-1)=real(CResistance(j-1));
      dUi(2*j-1,2*j)=-imag(CResistance(j-1));
      dUi(2*j,2*j-1)=imag(CResistance(j-1));
      dUi(2*j,2*j)=real(CResistance(j-1));
    }
    d U0=v.U(0);
    dUi(0,0)=zfd(mu0,rh);	// Zero frequency drag divided by zero-frequency velocity
    return dUi;
  }
  namespace laminardrag{
    // Resistance force for laminar flow for the zero-frequency. 
    d zerodrag_vert(d mu,d rh){
      return 3*mu/pow(rh,2);
    }
    d zerodrag_circ(d mu,d rh){
      return 2*mu/pow(rh,2);
    }
    d zerodrag_blapprox(d mu,d rh){ return 0; }
    
    ZeroFreqDrag::ZeroFreqDrag(const Tube& t){
      TRACE(0,"ZeroFreqDrag::ZeroFreqDrag()");
      if(t.geom.shape.compare("vert")==0)
	zerodrag_funptr=&zerodrag_vert;
      else if(t.geom.shape.compare("circ")==0)
	zerodrag_funptr=&zerodrag_circ;
      else if(t.geom.shape.compare("blapprox")==0)
	zerodrag_funptr=&zerodrag_blapprox;
      else
	{
	  WARN("Warning: tube.geom.shape unknown for ZeroFreqDrag. Aborting...");
	  abort();
	}
    }
    d ZeroFreqDrag::operator()(d mu,d rh,d U) const {
      return (*zerodrag_funptr)(mu,rh)*U;
    }
    d ZeroFreqDrag::operator()(d mu,d rh) const {
      return (*zerodrag_funptr)(mu,rh);
    }
  } // namespace laminardrag
} // namespace tube








