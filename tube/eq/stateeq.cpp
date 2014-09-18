#include "stateeq.h"
#include "tubevertex.h"

#define STATE_SCALE (1/v.gc->p0)
// #define STATE_SCALE (1)
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
  JacCol State::dpL(const TubeVertex& v) const {
    TRACE(0,"State::dpi");
    return JacCol(v.pL(),0.5*STATE_SCALE*eye<dmat>(v.gc->Ns,v.gc->Ns));
  }
  JacCol State::dpR(const TubeVertex& v) const {
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

  vd StateL::error(const TubeVertex& v)
   const {
    TRACE(6,"StateL::Error()");
    vd error(v.gc->Ns,fillwith::zeros);
    // TRACE(-1,"StateL p0:"<<p0);
    error+=v.pL()();
    error(0)+=v.gc->p0;	       // Add p0 part
    vd rhoTti=v.rho.tdata()%v.T.tdata();
    error+=v.gc->gas.Rs()*(v.gc->fDFT*(v.sLWi*rhoTti));
    if(v.left)  
      error+=v.gc->gas.Rs()*(v.gc->fDFT*(v.sLWim1*v.left->rho.tdata()%v.left->T.tdata()));
    else
      error+=v.gc->gas.Rs()*(v.gc->fDFT*(v.sLWip1*v.right->rho.tdata()%v.right->T.tdata()));      
    return STATE_SCALE*error;
  }
  JacRow StateL::jac(const TubeVertex& v) const{
    TRACE(6,"StateL::jac()");
    JacRow jac(dofnr,5);
    jac+=dpL(v);
    jac+=dTi(v);
    jac+=drhoi(v);
    if(v.left){
      jac+=dTim1(v);
      jac+=drhoim1(v);
    }
    else{
      jac+=dTip1(v);
      jac+=drhoip1(v);
    }
    return jac;
  }
  JacCol StateL::dpL(const TubeVertex& v)
   const {
    TRACE(0,"StateL::dpi");
    return JacCol(v.pL(),STATE_SCALE*eye<dmat>(v.gc->Ns,v.gc->Ns));
  }

  JacCol StateL::dTi(const TubeVertex& v)
   const {
    TRACE(0,"StateL::dTi()");
    dmat rhotidiag=diagmat(v.rho.tdata());
 return  JacCol(v.T,v.sLWi*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*rhotidiag*v.gc->iDFT);
  }
  JacCol StateL::drhoi(const TubeVertex& v)
   const {
    TRACE(0,"StateL::drhoi()");
    dmat Ttidiag=diagmat(v.T.tdata());
    return JacCol(v.rho,v.sLWi*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*Ttidiag*v.gc->iDFT);
  }
  JacCol StateL::dTim1(const TubeVertex& v)
   const {
    TRACE(0,"StateL::dTim1()");
    dmat rhotim1diag=diagmat(v.left->rho.tdata());
    return JacCol(v.left->T,v.sLWim1*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*rhotim1diag*v.gc->iDFT);
  }
  JacCol StateL::drhoim1(const TubeVertex& v)
   const {
    TRACE(0,"StateL::drhoim1()");
    dmat Ttim1diag=diagmat(v.left->T.tdata());
    return JacCol(v.left->rho,v.sLWim1*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*Ttim1diag*v.gc->iDFT);
  }
  JacCol StateL::dTip1(const TubeVertex& v)
   const {
    TRACE(0,"StateL::dTip1()");
    dmat rhotip1diag=diagmat(v.right->rho.tdata());
    return JacCol(v.right->T,v.sLWip1*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*rhotip1diag*v.gc->iDFT);
  }
  JacCol StateL::drhoip1(const TubeVertex& v)
   const {
    TRACE(0,"StateL::drhoip1()");
    dmat Ttip1diag=diagmat(v.right->T.tdata());
    return JacCol(v.right->rho,v.sLWip1*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*Ttip1diag*v.gc->iDFT);
  }
    
  // *********************************************** Now for pR,
  // *********************************************** generally a last
  // *********************************************** node

    vd StateR::error(const TubeVertex& v)
   const {
    TRACE(6,"StateR::Error()");
    vd error(v.gc->Ns,fillwith::zeros);
    // TRACE(-1,"StateL p0:"<<p0);
    error+=v.pR()();
    error(0)+=v.gc->p0;	       // Add p0 part
    vd rhoTti=v.rho.tdata()%v.T.tdata();
    error+=-v.gc->gas.Rs()*(v.gc->fDFT*(v.w.wRNm1*rhoTti));
    if(v.left)  
      error+=-v.gc->gas.Rs()*(v.gc->fDFT*(v.w.wRNm2*v.left->rho.tdata()%v.left->T.tdata()));
    else{
      WARN("StateR equation not implemented correctly");
      exit(1);
    }
    return STATE_SCALE*error;
  }
  JacRow StateR::jac(const TubeVertex& v) const{
    TRACE(6,"StateR::jac()");
    JacRow jac(dofnr,5);
    jac+=dpR(v);
    jac+=dTi(v);
    jac+=drhoi(v);
    if(v.left){
      jac+=dTim1(v);
      jac+=drhoim1(v);
    }
    else{
      WARN("StateR equation not implemented correctly");
      exit(1);
    }
    return jac;
  }
  JacCol StateR::dpR(const TubeVertex& v)
   const {
    TRACE(0,"StateR::dpi");
    return JacCol(v.pR(),STATE_SCALE*eye<dmat>(v.gc->Ns,v.gc->Ns));
  }
  JacCol StateR::dTi(const TubeVertex& v)
   const {
    TRACE(0,"StateR::dTi()");
    dmat rhotidiag=diagmat(v.rho.tdata());
    return JacCol(v.T,-v.w.wRNm1*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*rhotidiag*v.gc->iDFT);
  }    
  JacCol StateR::drhoi(const TubeVertex& v)
   const {
    TRACE(0,"StateR::drhoi()");
    dmat Ttidiag=diagmat(v.T.tdata());
    return JacCol(v.rho,-v.w.wRNm1*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*Ttidiag*v.gc->iDFT);
  }

  JacCol StateR::dTim1(const TubeVertex& v)
   const {
    TRACE(0,"StateR::dTim1()");
    dmat rhotim1diag=diagmat(v.left->rho.tdata());
    return JacCol(v.left->T,-v.w.wRNm2*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*rhotim1diag*v.gc->iDFT);
  }
  JacCol StateR::drhoim1(const TubeVertex& v)
   const {
    TRACE(0,"StateR::drhoim1()");
    dmat Ttim1diag=diagmat(v.left->T.tdata());
    return JacCol(v.left->rho,-v.w.wRNm2*STATE_SCALE*v.gc->gas.Rs()*v.gc->fDFT*Ttim1diag*v.gc->iDFT);
  }    
} // namespace tube




