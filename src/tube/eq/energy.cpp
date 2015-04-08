// energy.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////

// #define TRACERPLUS 20
// #define ENERGY_SCALE (1/v.gc->rho0/v.gc->c0)
// #define ENERGY_SCALE (1.0/v.gc->p0)
// #define ENERGY_SCALE (1.0/100)
#define ENERGY_SCALE (1.0)
// #define TRACERPLUS 15
#ifdef NOHEAT
#error Noheat already defined!
#endif

#include "energy.h"
#include "cell.h"
#include "weightfactors.h"
#include "tube.h"
#include "jacrow.h"
#include <tuple>
#define NOHEAT

#define iDFT (v.gc->iDFT)
#define fDFT (v.gc->fDFT)
#define DDTfd (v.gc->DDTfd)
#define Ns (v.gc->Ns())

namespace tube{
  using tasystem::Jacobian;
  using tasystem::JacRow;
  using tasystem::JacCol;

  void Energy::show() const{
    cout << "------------- Full energy equation\n";
    cout << "Wddt     :"<<Wddt<<"\n";
    cout << "WLl     :"<<WLl<<"\n";
    cout << "WLr     :"<<WLr<<"\n";
    cout << "WRl     :"<<WRl<<"\n";
    cout << "WRr     :"<<WRr<<"\n";
    // cout << "WkinLl   :"<<WkinLl<<"\n";
    // cout << "WkinLr   :"<<WkinLr<<"\n";
    // cout << "WkinRl   :"<< WkinRl<<"\n";
    // cout << "WkinRr   :"<< WkinRr<<"\n";
    cout << "WcLl     :"<<WcLl<<"\n";
    cout << "WcLr     :"<<WcLr<<"\n";
    cout << "WcRl     :"<<WcRl<<"\n";
    cout << "WcRr     :"<<WcRr<<"\n";

  }
  void Energy::init(){
    TRACE(8,"Energy::init(tube)");
    const Tube& t=v.getTube();
    heat=&t.getHeatSource();
    std::tie(WLl,WLr,WRl,WRr)=WeightFactors(v);
    d vx=v.vx;
    d xL=v.xL;
    d xR=v.xR;    

    Wddt=v.vVf;
    // VARTRACE(25,v.vVf);
    // Difference of factor half
    Wddtkin=0.5*(v.xR-v.xL);

    // WkinL=-0.5*w.wLl/SfLsq;
    // Wkin=0.5*(w.wRl/SfRsq-w.wLr/SfLsq);
    // WkinR=0.5*w.wRr/SfRsq;
   
    // vxm1=0 for the leftmost cell, so this is always true:
    // WcLl=w.SfL/(w.vx-w.vxm1);
    // WcLr=-WcLl;
    TRACE(8,"Energy::init(tube)");

    // WkinL=-0.5*w.wLl/SfLsq;
    // Wkin=0.5*(w.wRl/SfRsq-w.wLr/SfLsq);
    // WkinR=0.5*w.wRr/SfRsq;
    

    if(v.left()){
      d vxim1=v.left()->vx;
      WcLl=v.SfL/(vx-vxim1);
      //   d vSfLsq=pow(w.vSfL,2);
      //   WkinLl=0.5*w.wLl/vSfLsq;
      //   WkinLr=0.5*w.wLr/vSfsq;
    }
    else{
      //   d SfLsq=pow(w.SfL,2);
      WcLl=v.SfL/vx;
      WcLr=-WcLl;
      //   WkinLl=0.5/SfLsq;
      //   WkinLr=0;
    }
    WcLr=-WcLl;
    
    if(v.right()){    
      d vxip1=v.right()->vx;
      //   d vSfRsq=pow(w.vSfR,2);

      WcRl=v.SfR/(vxip1-vx);
      //   WkinRl=0.5*(w.wRl/vSfsq);
      //   WkinRr=0.5*(w.wRr/vSfRsq);
    }
    else{
    //   d SfRsq=pow(w.SfR,2);
      WcRl=v.SfR/(v.xR-vx);
    //   WkinRl=0;
    //   WkinRr=0.5/SfRsq;
    }
    WcRr=-WcRl;
  }
  d Energy::Htot() const{
    TRACE(10,"Energy::Htot()");
    // vd H=0.5*(HL()+HR());
    // return H(0);
    return 0;
  }
  vd Energy::error() const {		// Error in momentum equation
    TRACE(15,"Energy::Error(), i="<<v.geti());
    assert(v.gc!=nullptr);

    // Time derivative of static enthalpy in cell
    d gamma=this->gamma();
    vd error=(Wddt*DDTfd*v.p()())/(gamma-1.0);
    // Time derivative of total energy
    // error+=Wddtkin*DDTfd*v.mu()();
    // Total enthalpy flux
    error+=v.mHR()()-v.mHL()();

    error+=QR()-QL();

    // External heat    
    assert(heat!=nullptr);
    #ifndef NOHEAT
    // error+=Wddt*heat->heat(v);
    #else
    if(v.geti()==0)
      WARN("Applying no heat coupling");
    #endif

    // (Boundary source term)
    // error+=v.esource();
    return error;
  }
  JacRow Energy::jac() const{
    TRACE(6,"Energy::jac()");
    JacRow jac(dofnr,12);
    d gamma=this->gamma();
    // Time-derivative of internal energy
    jac+=JacCol(v.p(),(Wddt*DDTfd)/(gamma-1.0));
    // jac+=JacCol(v.mu(),Wddtkin*DDTfd);

    // Enthalpy flow out minus in
    jac+=JacCol(v.mHR(),eye());
    jac+=JacCol(v.mHL(),-eye());

    // Heat conduction part
    jac+=dQR();
    jac+=(dQL()*=-1);

    // Transverse heat transver
    #ifndef NOHEAT
    // jac+=JacCol(v.U(),Wddt*heat->dUi(v));
    // jac+=JacCol(v.T(),Wddt*heat->dTi(v));
    #endif

    // jac*=ENERGY_SCALE;
    return jac;    
  }

  vd Energy::QL() const{
    TRACE(4,"Energy::QL()");
    const vd& Tt=v.T().tdata();
    const vd& Ttl=v.TL().tdata();
    // VARTRACE(15,fDFT*(kappaLt()%(WcLl*Ttl+WcLr*Tt)));
    return fDFT*(kappaLt()%(WcLl*Ttl+WcLr*Tt));
  }
  JacRow Energy::dQL() const{
    TRACE(4,"Energy::dQL()");
    JacRow dQL(2);
    dQL+=JacCol(v.T(),fDFT*diagmat(WcLr*kappaLt())*iDFT);
    dQL+=JacCol(v.TL(),fDFT*diagmat(WcLl*kappaLt())*iDFT);
    return dQL;
  }
  vd Energy::QR() const{
    TRACE(4,"Energy::QR()");
    const vd& Tt=v.T().tdata();
    const vd& Ttr=v.TR().tdata();
    VARTRACE(15,fDFT*(kappaRt()%(WcRl*Tt+WcRr*Ttr)));
    return fDFT*(kappaRt()%(WcRl*Tt+WcRr*Ttr));
  }
  JacRow Energy::dQR() const{
    TRACE(4,"Energy::dQR()");
    JacRow dQR(2);
    dQR+=JacCol(v.T(),fDFT*diagmat(WcRl*kappaRt())*iDFT);
    dQR+=JacCol(v.TR(),fDFT*diagmat(WcRr*kappaRt())*iDFT);
    return dQR;
  }
  vd Energy::kappaRt()  const {		// Returns thermal conductivity time domain data
    TRACE(5,"Energy::kappaRt()");
    const vd& Tt=v.T().tdata();
    const vd& TtR=v.TR().tdata();
    // VARTRACE(25,v.gc->gas().kappa(WRr*TtR+WRl*Tt));
    return v.gc->gas().kappa(WRr*TtR+WRl*Tt);
  }
  vd Energy::kappaLt()  const {		// Returns thermal conductivity time domain data
    TRACE(5,"Energy::kappaRt()");
    const vd& Tt=v.T().tdata();
    const vd& TtL=v.TL().tdata();    
    // VARTRACE(25,v.gc->gas().kappa(WLl*TtL+WLr*Tt));
    return v.gc->gas().kappa(WLl*TtL+WLr*Tt);
  }
  void Energy::domg(vd& domg_) const {
    TRACE(0,"Energy::domg()");
    assert(v.gc!=nullptr);

    d gamma=this->gamma();
    d gamfac=gamma/(gamma-1.0);
    WARN("Is this correct??");
    vd domg_full=Wddt*DDTfd*v.p()()/(gamma-1.0)/v.gc->getomg(); // Static
    domg_full+=Wddtkin*DDTfd*v.mu()()/v.gc->getomg();                                                        // enthalpy
                                                             // term
    // domg_.subvec(dofnr+1,dofnr+2)=domg_full.subvec(1,2);
    domg_.subvec(dofnr,dofnr+Ns-1)=domg_full;
  }
  namespace{
    using std::make_tuple;
    std::tuple<d,d> weightfactors(const Cell& v){
      TRACE(15,"anonymous weightfactors for extrapolation of mH");
      assert((!v.left() && v.right()) || (v.left() && !v.right()));
      if(!v.left()){
        // Leftmost node
        d xL=0;
        d xR=v.xR;
        d xRR=v.right()->xR;

        // W0: weight factor for contribution of quantity at right
        // cell wall for something at the left cell wall

        // W1: weight factor for contribution of quantity at right
        // cell wall of neighbouring cell to what happens at the left
        // cell wall of this cell
        d W1=-xR/(xRR-xR);
        d W0=1-W1;
        // VARTRACE(40,W0);
        // VARTRACE(40,W1);
        return make_tuple(W0,W1);
      }
      else{
        d xR=v.xR;
        d xL=v.xL;
        d xLL=v.left()->xL;
        d WR2=(xR-xL)/(xLL-xL);
        d WR1=1-WR2;
        return make_tuple(WR1,WR2);
      }
    }
  }

  vd Energy::extrapolateEnthalpyFlow() const{
    TRACE(15,"Energy::extrapolateEnthalpyFlow()");
    // Can only be called for leftmost or rightmost node
    assert((!v.left() && v.right()) || (v.left() && !v.right()));
    if(!v.left()){
      d W0,W1; std::tie(W0,W1)=weightfactors(v);
      VARTRACE(40,W0);
      VARTRACE(40,W1);
      return W0*v.mHR()()+W1*v.right()->mHR()();
    }
    else{
      d WR1,WR2; std::tie(WR1,WR2)=weightfactors(v);
      VARTRACE(40,WR1);
      VARTRACE(40,WR2);
      return WR1*v.mHL()()+WR2*v.left()->mHL()();
    }
  }
  JacRow Energy::dExtrapolateEnthalpyFlow() const {
    TRACE(15,"Energy::dExtrapolateEnthalpyFlow()");
    JacRow jac(2);
    if(!v.left()){
      d W0,W1; std::tie(W0,W1)=weightfactors(v);
      VARTRACE(40,W0);
      VARTRACE(40,W1);
      jac+=JacCol(v.mHR(),W0*eye());
      jac+=JacCol(v.right()->mHR(),W1*eye());
    }
    else{
      d WR1,WR2; std::tie(WR1,WR2)=weightfactors(v);
      jac+=JacCol(v.mHL(),WR1*eye());
      jac+=JacCol(v.left()->mHL(),WR2*eye());
    }
    return jac;
  }
  vd Energy::extrapolateHeatFlow() const {
    TRACE(5,"Energy::extrapolateHeatFlow()");
    vd Qb(Ns);
    if(!v.left()){
      vd kappaLt=this->kappaLt();
      Qb=fDFT*(kappaLt%(WcLl*v.TL().tdata()+WcLr*v.T().tdata()));
    }
    else if(!v.right()){
      vd kappaRt=this->kappaRt();
      Qb=fDFT*(kappaRt%(WcRl*v.T().tdata()+WcRr*v.TR().tdata()));
    }
    else{
      WARN("Extrapolation of heat flow asked for an inner cell. This is a bug. Aborting...");
      abort();
    }
    return Qb;
  }
  JacRow Energy::dExtrapolateHeatFlow() const{
    TRACE(5,"Energy::dExtrapolateHeatFlow()");
    JacRow dQb(2);

    // VARTRACE(30,kappaLt)      ;
    // VARTRACE(30,kappaRt);
    if(!v.left()){
      vd kappaLt=this->kappaLt();
      dQb+=JacCol(v.T(),fDFT*diagmat(WcLr*kappaLt)*iDFT);
      dQb+=JacCol(v.TL(),fDFT*diagmat(WcLl*kappaLt)*iDFT);
    }
    else if(!v.right()){
      vd kappaRt=this->kappaRt();
      dQb+=JacCol(v.T(),fDFT*diagmat(WcRl*kappaRt)*iDFT);
      dQb+=JacCol(v.TR(),fDFT*diagmat(WcRr*kappaRt)*iDFT);
    }
    else{
      WARN("That went fatally wrong!");
      abort();
    }
    return dQb;
  }
  d Energy::gamma() const {
    // d T0=v.T(0);
    d T0=v.gc->T0(); 
    return v.gc->gas().gamma(T0);
  }


} // namespace tube

//////////////////////////////////////////////////////////////////////

