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
#include "jacrow.h"
#include "jacobian.h"

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
  SimpleTubeConnector::SimpleTubeConnector(const SimpleTubeConnector& other,const TaSystem& sys):
    Connector(other,sys),
    segnrs(other.segnrs),
    pos(other.pos)
  {
    bccells[0]=&sys.getTube(segnrs[0]).bcCell(pos[0]);
    bccells[1]=&sys.getTube(segnrs[1]).bcCell(pos[1]);
    setInit(true);
  }

  vd SimpleTubeConnector::error() const{
    TRACE(10,"SimpleTubeConnector::error()");
    d out0=1,out1=1;
    if(pos[0]==Pos::left)
      out0=-1;
    if(pos[1]==Pos::left)
      out1=-1;

    vd error(Neq*Ns);
    {
      vd errorm(Ns,fillwith::zeros);
      vd errormint(Ns,fillwith::zeros);    

      errorm=out0*bccells[0]->mbc()();
      errorm+=out1*bccells[1]->mbc()();
      // Interpolation: mass flow is average of extrapolated mf
      errormint+=0.5*out0*bccells[0]->extrapolateQuant(massFlow);
      errormint+=0.5*out1*bccells[1]->extrapolateQuant(massFlow);
      errormint-=out0*bccells[0]->mbc()();

      error.subvec(0,Ns-1)=errorm;
      error.subvec(Ns,2*Ns-1)=errormint;
    }
    {
      
      vd errormH(Ns,fillwith::zeros);
      vd errormHint(Ns,fillwith::zeros);    

      errormH+=out0*bccells[0]->mHbc()();
      errormH+=out1*bccells[1]->mHbc()();
      
      errormHint+=0.5*out0*bccells[0]->extrapolateQuant(enthalpyFlow);
      errormHint+=0.5*out1*bccells[1]->extrapolateQuant(enthalpyFlow);
      errormHint-=out0*bccells[0]->mHbc()();
      
      error.subvec(2,3*Ns-1)=errormH;
      error.subvec(3*Ns,4*Ns-1)=errormHint;

    }
    {
      
      vd errorQ=out0*bccells[0]->extrapolateQuant(heatFlow);
      errorQ+=out1*bccells[1]->extrapolateQuant(heatFlow);

      // d Sfgem=0.5*(bccells[0]->Sfbc()+bccells[1]->Sfbc());
      vd errorM=(bccells[0]->extrapolateQuant(pressure)-
                 bccells[1]->extrapolateQuant(pressure));

      error.subvec(4,5*Ns-1)=errorQ;
      error.subvec(5*Ns,6*Ns-1)=errorM;

    }    

    return error;
  }
  void SimpleTubeConnector::jac(tasystem::Jacobian& jac) const {
    TRACE(15,"SimpleTubeconnector::jac()");
    d out0=1,out1=1;
    us eqnr=firsteqnr;
    
    if(pos[0]==Pos::left)
      out0=-1;
    if(pos[1]==Pos::left)
      out1=-1;
    {
      JacRow mjac(eqnr,2);
      eqnr+=Ns;
      mjac+=JacCol(bccells[0]->mbc(),out0*eye);
      mjac+=JacCol(bccells[1]->mbc(),out1*eye);
      
      jac+=mjac;
      JacRow mjacint(eqnr,5);
      eqnr+=Ns;
      // // Interpolation: mass flow is average of extrapolated mf
      mjacint+=(bccells[0]->dExtrapolateQuant(massFlow)*=0.5*out0);
      mjacint+=(bccells[1]->dExtrapolateQuant(massFlow)*=0.5*out1);
      mjacint+=JacCol(bccells[0]->mbc(),-out0*eye);

      jac+=mjacint;
    }
    {
      JacRow mHjac(eqnr,2);
      eqnr+=Ns;
      mHjac+=JacCol(bccells[0]->mHbc(),out0*eye);
      mHjac+=JacCol(bccells[1]->mHbc(),out1*eye);
      
      jac+=mHjac;
      JacRow mHjacint(eqnr,5);
      eqnr+=Ns;

      // // Interpolation: mass flow is average of extrapolated mf
      mHjacint+=(bccells[0]->dExtrapolateQuant(enthalpyFlow)*=0.5*out0);
      mHjacint+=(bccells[1]->dExtrapolateQuant(enthalpyFlow)*=0.5*out1);
      mHjacint+=JacCol(bccells[0]->mHbc(),-out0*eye);

      jac+=mHjacint;
    }    
    {
      JacRow Qjac(eqnr,5);
      eqnr+=Ns;

      Qjac+=(bccells[0]->dExtrapolateQuant(heatFlow)*=out0);
      Qjac+=(bccells[1]->dExtrapolateQuant(heatFlow)*=out1);
      jac+=Qjac;

      JacRow pjac(eqnr,2);
      // eqnr+=Ns; // Not needed no row below
      pjac+=bccells[0]->dExtrapolateQuant(pressure);
      pjac+=(bccells[1]->dExtrapolateQuant(pressure)*=-1);
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



