#include "stateeq.h"
#include "tubevertex.h"

#define STATE_SCALE (1.0/v.gc->p0)

namespace tube{

  vd State::error(const TubeVertex& v)
   const {
    TRACE(6,"State::Error()");
    vd error(v.gc->Ns,fillwith::zeros);
    // TRACE(-1,"State p0:"<<p0);
    error+=0.5*(v.pL()()+v.pR()());
    error(0)+=v.gc->p0;	       // Add p0 part
    // TRACE(-1,"state error:"<<error);    
    // TRACE(-1,"T0:"<<vertex.gc->gas.Rs()*fDFT()*(vertex.T.tdata()%vertex.rho.tdata()));    
    error+=-1.0*v.gc->gas.Rs()*v.gc->fDFT*(v.rho.tdata()%v.T.tdata());
    // TRACE(-1,"state error:"<<error);
    return STATE_SCALE*error;
  }
  JacRow State::jac(const TubeVertex& v) const{
    TRACE(6,"State::jac()");
    JacRow jac(dofnr,3);
    jac+=dpL(v);
    jac+=dpR(v);    
    jac+=dTi(v);
    jac+=drhoi(v);
    return jac;
  }
  JacCol State::dpL(const TubeVertex& v)
   const {
    TRACE(0,"State::dpi");
    return JacCol(v.pL(),0.5*STATE_SCALE*eye<dmat>(v.gc->Ns,v.gc->Ns));
  }
  JacCol State::dpR(const TubeVertex& v)
   const {
    TRACE(0,"State::dpi");
    return JacCol(v.pR(),0.5*STATE_SCALE*eye<dmat>(v.gc->Ns,v.gc->Ns));
  }
  JacCol State::dTi(const TubeVertex& v)
   const {
    TRACE(0,"State::dTi()");
    dmat rhotidiag=diagmat(v.rho.tdata());
    return JacCol(v.T,-1.0*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*rhotidiag*v.gc->iDFT);
  }
  JacCol State::drhoi(const TubeVertex& v)
   const {
    TRACE(0,"State::drhoi()");
    dmat Ttidiag=diagmat(v.T.tdata());
    return JacCol(v.rho,-1.0*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*Ttidiag*v.gc->iDFT);
  }

} // namespace tube




