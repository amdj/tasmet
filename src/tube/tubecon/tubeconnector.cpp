// tubeconnector.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
// Connect two tube-type segments by using conservation of mass,
// momentum and energy. Number of equations 
//////////////////////////////////////////////////////////////////////

#include "tubeconnector.h"
#include "tasystem.h"
#include "tube.h"
#include "bccell.h"
#include "constants.h"
#include "jacobian.h"
#include "energy.h"

#define fDFT (gc->fDFT)
#define iDFT (gc->iDFT)
#define DDTfd (gc->DDTfd)

#define Ns (gc->Ns())
#define eye (arma::eye(Ns,Ns))

namespace tube {
  using tasystem::TaSystem;
  using tasystem::JacRow;
  using tasystem::JacCol;
  // Number of equations corresponding to this connection:
  const int Neq=6;

  SimpleTubeConnector::SimpleTubeConnector(us seg1,Pos pos1,\
                                           us seg2,Pos pos2)
    
  {
    TRACE(15,"SimpleTubeConnector::SimpleTubeConnector()");
  
    if(max(seg1,seg2)>constants::maxsegs)
      throw MyError("Too high segment number given");
    if(seg1==seg2)
      throw MyError("Segments cannot be the same");
    segnrs[0]=seg1;
    segnrs[1]=seg2;
    pos[0]=pos1;
    pos[1]=pos2;
      
  }
  SimpleTubeConnector::SimpleTubeConnector(const SimpleTubeConnector& o,
                                           const TaSystem& sys):
    Connector(o,sys),
    segnrs(o.segnrs),
    pos(o.pos)
  {
    bccells[0]=&sys.getTube(segnrs[0]).bcCell(pos[0]);
    bccells[1]=&sys.getTube(segnrs[1]).bcCell(pos[1]);

    if(pos[0]==Pos::left) {
      out[0]=-1;
      dx+=bccells[0]->vx;      
    }
    else{
      dx+=bccells[0]->xR-bccells[0]->vx;
    }
    if(pos[1]==Pos::left) {
      out[1]=-1;
      dx+=bccells[1]->vx;      
    }
    else
      dx+=bccells[1]->xR-bccells[1]->vx;

    Sfgem=0.5*(bccells[0]->Sfbc()+bccells[1]->Sfbc());
    
    VARTRACE(15,pos[0]);
    VARTRACE(15,pos[1]);
    VARTRACE(15,out[0]);
    VARTRACE(15,out[1]);
    setInit(true);
  }
  vd SimpleTubeConnector::kappaSft() const  {
    TRACE(10,"SimpleTubeConnector::kappaSft()");
    d dx=0;
    // Distance from left interior node to right interior node
    vd kappaSf0t,kappaSf1t;

    if(pos[0]==Pos::left){
      kappaSf0t=Energy::kappaLt(*bccells[0]);

    }
    else{
      kappaSf0t=Energy::kappaRt(*bccells[0]);

    }
    if(pos[1]==Pos::left){
      kappaSf1t=Energy::kappaLt(*bccells[1]);
    }
    else{
      kappaSf1t=Energy::kappaRt(*bccells[1]);
    }
    return 0.5*(kappaSf0t*bccells[0]->Sfbc()
                  +kappaSf1t*bccells[1]->Sfbc());
  }


  vd SimpleTubeConnector::error() const{
    TRACE(10,"SimpleTubeConnector::error()");

    us nr=0;
    vd error(Neq*Ns);
    const vd& T1t=bccells[0]->T().tdata();
    const vd& T2t=bccells[1]->T().tdata();
      d T0=gc->T0();
      d cp=gc->gas().cp(T0);
      const vd p1t=bccells[0]->extrapolateQuant(Varnr::p);
      const vd p2t=bccells[2]->extrapolateQuant(Varnr::p);

    // First equation: mass balance
    error.subvec(Ns*nr,Ns*(nr+1)-1)=\
      out[0]*bccells[0]->mbc()()\
      +out[1]*bccells[1]->mbc()();
    nr++;
    // Second equation: enthalpy flow balance
    error.subvec(Ns*nr,Ns*(nr+1)-1)=\
      out[0]*bccells[0]->mHbc()()\
      +out[1]*bccells[1]->mHbc()();

    nr++;
      // Third equation: Enthalpy flow at the interface equals the average of the
      // enthalpy flows of both tubes
      error.subvec(Ns*nr,Ns*(nr+1)-1)=\
        0.5*out[0]*bccells[0]->extrapolateQuant(Varnr::mH)\
        -0.5*out[1]*bccells[1]->extrapolateQuant(Varnr::mH)\
        -out[0]*bccells[0]->mHbc()();

      nr++;
    {
      vd kappaSft=this->kappaSft();
      
      error.subvec(Ns*nr,Ns*(nr+1)-1)=\
        fDFT*(kappaSft%(T1t-T2t)/dx)\
        -out[0]*bccells[0]->extrapolateQuant(Varnr::Q);
      nr++;      
    }
    {
      error.subvec(Ns*nr,Ns*(nr+1)-1)=\
        out[0]*bccells[0]->extrapolateQuant(Varnr::Q)\
        +out[1]*bccells[1]->extrapolateQuant(Varnr::Q);

      nr++;
    }
    {
      // Change in exergy
      // DeltaE=T0*cp*log(
      error.subvec(Ns*nr,Ns*(nr+1)-1)=\
        Sfgem*(bccells[1]->extrapolateQuant(Varnr::p)\
               -bccells[0]->extrapolateQuant(Varnr::p));//+
        // bccells[1]->extrapolateQuant(Varnr::mu)
        // -bccells[0]->extrapolateQuant(Varnr::mu);
      nr++;
    }    

    return error;
  }
  void SimpleTubeConnector::jac(tasystem::Jacobian& jac) const {
    TRACE(15,"SimpleTubeconnector::jac()");
    us eqnr=firsteqnr;
    
    {                           // Mass flow continuity
      JacRow mjac(eqnr,2);
      eqnr+=Ns;
      mjac+=JacCol(bccells[0]->mbc(),out[0]*eye);
      mjac+=JacCol(bccells[1]->mbc(),out[1]*eye);
      jac+=mjac;
    }
    {
      JacRow mHjac(eqnr,2);
      eqnr+=Ns;
      mHjac+=JacCol(bccells[0]->mHbc(),out[0]*eye);
      mHjac+=JacCol(bccells[1]->mHbc(),out[1]*eye);
      
      jac+=mHjac;
    }
    {
      JacRow mHjacint(eqnr,5);
      eqnr+=Ns;

      mHjacint+=(bccells[0]->dExtrapolateQuant(Varnr::mH)*=0.5*out[0]);
      mHjacint+=(bccells[1]->dExtrapolateQuant(Varnr::mH)*=-0.5*out[1]);
      mHjacint+=JacCol(bccells[0]->mHbc(),-out[0]*eye);

      jac+=mHjacint;
    }
    {
      JacRow Qjac(eqnr,5);
      eqnr+=Ns;

      vd kappaSft=this->kappaSft();
      VARTRACE(15,kappaSft);
      Qjac+=JacCol(bccells[0]->T(),fDFT*diagmat(kappaSft/dx)*iDFT);
      Qjac+=JacCol(bccells[1]->T(),-fDFT*diagmat(kappaSft/dx)*iDFT);
      Qjac+=(bccells[0]->dExtrapolateQuant(Varnr::Q)*=-out[0]);

      jac+=Qjac;

    }    
    {
      JacRow Qintjac(eqnr,5);
      eqnr+=Ns;

      Qintjac+=(bccells[0]->dExtrapolateQuant(Varnr::Q)*=out[0]);
      Qintjac+=(bccells[1]->dExtrapolateQuant(Varnr::Q)*=out[1]);
      jac+=Qintjac;
    }
    {
      JacRow pjac(eqnr,2);
      // eqnr+=Ns; // Not needed no row below
      pjac+=(bccells[1]->dExtrapolateQuant(Varnr::p)*=Sfgem);
      pjac+=(bccells[0]->dExtrapolateQuant(Varnr::p)*=-Sfgem);
      // pjac+=(bccells[1]->dExtrapolateQuant(Varnr::mu));
      // pjac+=(bccells[0]->dExtrapolateQuant(Varnr::mu)*=-1);             
      jac+=pjac;
    }    

  }
  void SimpleTubeConnector::setEqNrs(us firsteqnr){
    TRACE(15,"SimpleTubeConnector::setEqNrs()");
    this->firsteqnr=firsteqnr;
  }
  us SimpleTubeConnector::getNEqs() const {
    TRACE(15,"SimpleTubeConnector::getNEqs()");
    return Neq*Ns;
  }
  void SimpleTubeConnector::show(us) const{
    TRACE(15,"SimpleTubeConnector::show()");
    checkInit();
    
    cout << "SimpleTubeConnector which connects tube " << segnrs[0] <<
      " at the " << posWord(pos[0]) << " side to tube " << segnrs[1] <<
      " on the " << posWord(pos[1]) << " side.\n";
  }
  void SimpleTubeConnector::updateNf(){
    TRACE(15,"SimpleTubeConnector::updateNf()");
  }

} // namespace tube
//////////////////////////////////////////////////////////////////////



