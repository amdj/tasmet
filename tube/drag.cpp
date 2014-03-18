#include "drag.h"
#include "tube.h"


namespace tube{

  DragResistance::DragResistance(const Tube& t):tube(t),
vop(tube.vop)  {  }
  dmat DragResistance::dUi(us i) const {
    const us& Ns=vop.Ns;
    return dmat(Ns,Ns,fillwith::zeros); 
  }
  DragResistance::~DragResistance(){}
  vd  DragResistance::operator()(us i) const {
    //Compute the drag resistance for each frequency at node i in the momentum equation
    vd D(vop.Ns,fillwith::zeros);
    return D;
  }
vc LaminarDragResistance::ComplexResistancecoef(us i) const {
    TRACE(0,"LaminarDragResistance::ComplexResistancecoef()");
    const us& Nf=vop.Nf;
    d T0=tube.vvertex[i].T(0);	// Time-averaged temperature
    d mu0=tube.gas.mu(T0);
    d p0=tube.vvertex[i].p(0);
    d rho0=tube.gas.rho(T0,p0);

    const d& rh=tube.geom.rh(i);
    TRACE(-1,"rh: "<<rh);
    d omg,deltanu; 
    vc rh_over_deltanu(Nf);
    vc omgvec(Nf);
    for(us j=0;j<vop.Nf;j++)
      {
	TRACE(-1,"j:"<<j);
	omg=vop.omg*(j+1);
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
  LaminarDragResistance::LaminarDragResistance(const Tube&t):DragResistance(t),zfd(t){
    rf=math_common::rottfuncs(tube.geom.shape); // Reinitialize thermoviscous functions with right shape
  }
  vd LaminarDragResistance::operator()(us i) const {

    const d& rh=tube.geom.rh(i);
    TRACE(0,"LaminarDragResistance::operator()");
    const us& Nf=vop.Nf;
    dmat dDdU=dUi(i);
    vd drag=dDdU*tube.vvertex[i].U();
    return drag;
  }
  dmat LaminarDragResistance::dUi(us i) const { // Derivative of drag resistance to velocity
    TRACE(0,"LaminarDragResistance::dUi()");
    dmat dUi(vop.Ns,vop.Ns,fillwith::zeros);
    const d& rh=tube.geom.rh(i);
    const us& Nf=vop.Nf;
    d T0=tube.vvertex[i].T(0);	// Time-averaged temperature
    d mu0=tube.gas.mu(T0);

    // The complex resistance coefficient vector has size Nf
    vc CResistance=ComplexResistancecoef(i);
    for(us j=1;j<Nf+1;j++){
      dUi(2*j-1,2*j-1)=real(CResistance(j-1));
      dUi(2*j-1,2*j)=-imag(CResistance(j-1));
      dUi(2*j,2*j-1)=imag(CResistance(j-1));
      dUi(2*j,2*j)=real(CResistance(j-1));
    }
    d U0=tube.vvertex[i].U(0);
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
    ZerofreqDrag::ZerofreqDrag(const Tube& t): tube(t){
      TRACE(0,"ZerofreqDrag::ZerofreqDrag()");
      if(tube.geom.shape=="vert")
	zerodrag_funptr=&zerodrag_vert;
      else if(tube.geom.shape=="circ")
	zerodrag_funptr=&zerodrag_circ;
      else if(tube.geom.shape=="blapprox")
	zerodrag_funptr=&zerodrag_blapprox;
      else
	{
	  TRACE(10,"Warning: tube.geom.shape unknown for zerofreqdrag, exiting");
	  exit(1);
	}
    }
    ZerofreqDrag::~ZerofreqDrag(){}
    d ZerofreqDrag::operator()(d mu,d rh,d U) const {
      return (*zerodrag_funptr)(mu,rh)*U;
    }
    d ZerofreqDrag::operator()(d mu,d rh) const {
      return (*zerodrag_funptr)(mu,rh);
    }
  } // namespace laminardrag
} // namespace tube
