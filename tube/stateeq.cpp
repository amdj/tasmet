#include "stateeq.h"
#include "vertex.h"
#include "tube.h"
#include "../var/var.h"
namespace tube{

  State::State(const Tube& tube,TubeVertex& gp):Equation(tube,gp){
    TRACE(0,"State constructor done");
  }

  vd State::Error()
  {
    TRACE(0,"State::Error()");
    vd error(Ns,fillwith::zeros);
    vd p0=getp0();
    // TRACE(-1,"State p0:"<<p0);
    error+=p0+vertex.p();
    // TRACE(-1,"state error:"<<error);    
    // TRACE(-1,"T0:"<<tube.gas.Rs()*fDFT*(vertex.T.tdata()%vertex.rho.tdata()));    
    error+=-1.0*tube.gas.Rs()*fDFT*(vertex.rho.tdata()%vertex.T.tdata());
    // TRACE(-1,"state error:"<<error);
    return error;
  }
  dmat State::dpi()
  {
    TRACE(0,"State::dpi");
    return eye<dmat>(Ns,Ns);
  }
  dmat State::dTi()
  {
    TRACE(0,"State::dTi()");
    dmat rhotidiag=diagmat(vertex.rho.tdata());
    return -1.0*tube.gas.Rs()*vop.fDFT*rhotidiag*vop.iDFT;
  }
  dmat State::drhoi()
  {
    TRACE(0,"State::drhoi()");
    dmat Ttidiag=diagmat(vertex.T.tdata());
    return -1.0*tube.gas.Rs()*vop.fDFT*Ttidiag*vop.iDFT;
  }
  State::~State(){
    TRACE(-5,"State eq destructor");
  }
} // namespace tube
