// file: bcvertex.cpp, created March 20th, 2014.
// Author: J.A. de Jong
#include "bcvertex.h"
#include "momscale.h"


namespace tube{

  TubeBcVertex::TubeBcVertex(const Tube& t,us vertexnr):TubeVertex(t,vertexnr){

  }  
  TubeBcVertex::~TubeBcVertex(){}
  LeftPressure::LeftPressure(const Tube& t):TubeBcVertex(t,0),pL(gc),TL(gc){
    pL(0)=tube.gc.p0;
    TL(0)=tube.gc.T0;
    Init();
   }
   LeftPressure::LeftPressure(const Tube& t,variable::var& pres,variable::var& temp):TubeBcVertex(t,0),pL(pres),TL(temp) {
     TRACE(0,"LeftPressure full constructor");
     Init();
   }
   LeftPressure::LeftPressure(const Tube& t,variable::var& pres):TubeBcVertex(t,0),pL(pres),TL(gc){
     TRACE(0,"LeftPressure constructor for given pressure. Temperature computed");    
     d T0=tube.gc.T0;
     d gamma=tube.gas.gamma(T0);
     vd p0(Ns,fillwith::ones); p0*=tube.gc.p0;
     // TRACE(-1,"p0:"<<p0);
     vd TLt=T0*pow((p0+pL.tdata())/p0,(gamma-1.0)/gamma);		// Adiabatic compression/expansion
     // TRACE(-1,"TLt:"<<TLt);
     TL.settdata(TLt);

     // TRACE(-1,"rho:"<<rho());
     Init();
   }
   void LeftPressure::Init(){
     // Change continuity equation for an open boundary
     c.Wi=c.wRl-c.wL0;
     c.Wip1=c.wRr-c.wL1;

     // Change momentum equation for open boundary, and prescribed pressure
     m.Wuim1=0;
     m.Wui=m.wRl-m.wL0;
     m.Wuip1=m.wRr-m.wL1;
     m.Wpi=m.wRl*m.SfR+m.SfL-m.SfR;
     m.Wpip1=m.SfR*m.wRr;
     // Change energy equation for open boundary and prescribed pressure
     e.Wgim1=0;
     e.Wgi=e.wRl-e.wL0;
     e.Wgip1=e.wRr-e.wL1;

     e.Wjim1=0;
     e.Wji=e.wL0-e.wRl;
     e.Wjip1=e.wL1-e.wRr;

     d vxi=tube.geom.vx(0);
     d vxip1=tube.geom.vx(1);
     d dxp=vxip1-vxi;

     e.Wc1=0;
     e.Wc2=e.SfL/vxi;
     e.Wc3=e.SfR/dxp;
     e.Wc4=-e.SfR/dxp;
       // TODO Fill this further!

   }
   vd LeftPressure::msource() const{
     TRACE(0,"LeftPressure::msource()");
     vd msource(Ns,fillwith::ones);
     msource=-1.0*e.SfL*pL();
 // This one should not yet be scaled. The scaling is done in the
 // error term after adding this source.
    TRACE(-1,"msource:"<<msource);
    return msource;
  }
  vd LeftPressure::esource() const {
    TRACE(0,"LeftPressure::esource()");
    vd esource(Ns,fillwith::zeros);
    d vxi=tube.geom.vx(0);
    vd TLt=TL.tdata();
    vd kappaL=tube.gas.kappa(TLt);
    // TRACE(10,"Important: put esource on when going back to full energy eq!");
    esource+=-1.0*e.SfL*(e.fDFT*(kappaL%TLt))/vxi;
    TRACE(-1,"esource:"<<esource);
    return esource;    
  }
  LeftPressure::~LeftPressure(){
    TRACE(-5,"LeftPressure destructor");
  }

  RightImpedanceMomentumEq::RightImpedanceMomentumEq(const Tube& t,TubeBcVertex& tv,vd& Z):Momentum(t,tv),Z(Z){
    TRACE(0,"RightImpedanceMomentumEq::RightImpedanceMomentumEq()");
  }
  vd RightImpedanceMomentumEq::Error(){
    TRACE(0,"RightImpedanceMomentumEq::Error()");
    vd error(Ns,fillwith::zeros);
    // Add the normal stuff
    error+=Momentum::Error();  
    
    // And add the right pressure as being an impedance times the
    // velocity:
    vd errorZ=MOM_SCALE*SfR*Z%(wRNm1*vertex.U()+wRNm2*left->U());

    // errorZ(0)*=MOM_SCALE0;
    error+=errorZ;


    // Use time-averaged momentum equation to fix the pressure to p0
    // at the last node
    // Impedance-like
    // error(0)=1.0*MOM_SCALE0*MOM_SCALE*(vertex.p(0)-Z(0)*vertex.U(0));

    // Just set velocity at last node equalt to zero :Velocity zero
    // error(0)=1.0*MOM_SCALE0*MOM_SCALE*vertex.U(0);

    // Pressure opening
    // error(0)+=wRNm1*MOM_SCALE0SfR*vertex.p(0)+wRNm2*SfR*left->p(0);

    
    // Pressure zero
    error(0)=MOM_SCALE0*MOM_SCALE*vertex.p(0);


    return error;
  }
  dmat RightImpedanceMomentumEq::dpi(){
    dmat dpi=Momentum::dpi();
    dpi.row(0).zeros();

    // Set mean pressure to zero on last node
    dpi(0,0)=MOM_SCALE*MOM_SCALE0;
    return dpi;
  }
  dmat RightImpedanceMomentumEq::dUi(){
    TRACE(0,"RightImpedanceMomentumEq::dUi()");
    dmat dUi=Momentum::dUi();
    dmat adddUi=MOM_SCALE*wRNm1*SfR*diagmat(Z);
    adddUi.row(0)*=MOM_SCALE0;
    dUi+=adddUi;

    // For pressure boundary condition
    dUi.row(0).zeros();
    
    // For velocity boundary condition
    // dUi.row(0).zeros();
    // dUi(0,0)=MOM_SCALE0*MOM_SCALE;
      
    return dUi;
  }

  dmat RightImpedanceMomentumEq::dUim1(){
    TRACE(0,"RightImpedanceMomentumEq::dUim1()");
    dmat dUim1=Momentum::dUim1();    
    dmat adddUim1=MOM_SCALE*wRNm2*SfR*diagmat(Z);
    adddUim1.row(0)*=MOM_SCALE0;
    dUim1+=adddUim1;

    // For velocity boundary condition
    dUim1.row(0).zeros();
    
    return dUim1;
  }

  dmat RightImpedanceMomentumEq::dpim1(){
    dmat dpim1=Momentum::dpim1();

    // For velocity and pressure boundary condition
    dpim1.row(0).zeros();

    return dpim1;
  }
  dmat RightImpedanceMomentumEq::drhoi(){
    TRACE(0,"RightImpedanceMomentumEq::drhoi()");
    dmat drhoi=Momentum::drhoi();

    // For velocity and pressure boundary condition
    drhoi.row(0).zeros();
    return drhoi;
  }

  dmat RightImpedanceMomentumEq::drhoim1(){
    TRACE(0,"RightImpedanceMomentumEq::drhoim1()");
    dmat drhoim1=Momentum::drhoim1();

    // For velocity and pressure boundary condition
    drhoim1.row(0).zeros();
    return drhoim1;
  }

  
  RightImpedanceMomentumEq::~RightImpedanceMomentumEq(){}
  RightImpedance::RightImpedance(const Tube& t,vd Z1):TubeBcVertex(t,t.Ncells-1),Z(Z1),mright(t,*this,Z){
    TRACE(-5,"RightImpedance constructor");
    // Change continuity equation for open boundary
    c.Wim1=c.wRNm2-c.wLl;
    c.Wi  =c.wRNm1-c.wLr;
    c.Wip1=0;

    // Change momentum eq for open boundary
    mright.Wuim1=-c.wLl/c.SfL+c.wRNm2/c.SfR;
    mright.Wui	=-c.wLr/c.SfL+c.wRNm1/c.SfR;
    mright.Wuip1=0;

    mright.Wpim1=-c.SfL*c.wLl;
    mright.Wpi	=-c.SfL*c.wLr+(m.SfL-m.SfR);
    mright.Wpip1=0;
    
    e.Wgim1=-e.wLl+e.wRNm2;
    e.Wgi  =-e.wLr+e.wRNm1;
    e.Wgip1=0;
    
    e.Wjim1=e.wLl-e.wRNm2;
    e.Wji  =e.wLr-e.wRNm1;
    e.Wjip1=0;
    
    d xi=tube.geom.vx(i);
    d xim1=tube.geom.vx(i-1);	 
    d dxm=xi-xim1;
    e.Wc1=-e.SfL/dxm;
    e.Wc2= e.SfL/dxm;
    e.Wc3=0;
    e.Wc4=0;

    // Last but not least: point momentum eq to new equation!
    eq[1]=&mright;
    // Conduction terms are not changed.
    
  }
  RightImpedance::~RightImpedance(){
    TRACE(-5,"RightImpedance destructor");
  }

} // namespace tube












