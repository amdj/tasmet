#include "stateeq.h"
#include "tubevertex.h"

#define STATE_SCALE (1.0)

namespace tube{

  State::State(TubeVertex& gp):TubeEquation(gp){
    TRACE(0,"State constructor done");
  }

  vd State::Error()
  {
    TRACE(6,"State::Error()");
    vd error(vertex.gc->Ns,fillwith::zeros);
    vd p0=getp0();
    // TRACE(-1,"State p0:"<<p0);
    error+=(p0+vertex.p());
    // TRACE(-1,"state error:"<<error);    
    // TRACE(-1,"T0:"<<vertex.gc->gas.Rs()*fDFT()*(vertex.T.tdata()%vertex.rho.tdata()));    
    error+=-1.0*vertex.gc->gas.Rs()*vertex.gc->fDFT*(vertex.rho.tdata()%vertex.T.tdata());
    // TRACE(-1,"state error:"<<error);
    return STATE_SCALE*error;
  }
  dmat State::dpi()
  {
    TRACE(0,"State::dpi");
    return STATE_SCALE*eye<dmat>(vertex.gc->Ns,vertex.gc->Ns);
  }
  dmat State::dTi()
  {
    TRACE(0,"State::dTi()");
    dmat rhotidiag=diagmat(vertex.rho.tdata());
    return -1.0*STATE_SCALE*vertex.gc->gas.Rs()*vertex.gc->fDFT*rhotidiag*vertex.gc->iDFT;
  }
  dmat State::drhoi()
  {
    TRACE(0,"State::drhoi()");
    dmat Ttidiag=diagmat(vertex.T.tdata());
    return -1.0*STATE_SCALE*vertex.gc->gas.Rs()*vertex.gc->fDFT*Ttidiag*vertex.gc->iDFT;
  }
  State::~State(){
    TRACE(-5,"State eq destructor");
  }
} // namespace tube
