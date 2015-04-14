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
    us nr=0;
    vd error(Neq*Ns);
    {
      vd errorm=out0*bccells[0]->mbc()();
      +out1*bccells[1]->mbc()();
      error.subvec(Ns*nr,Ns*(nr+1)-1)=errorm; nr++;
    }
    {
      vd errormH=out0*bccells[0]->mHbc()()
        +out1*bccells[1]->mHbc()();

      error.subvec(Ns*nr,Ns*(nr+1)-1)=errormH;  nr++;
    }
    {
      // Interpolation: mass flow is average of extrapolated mf
      vd errormHint=0.5*out0*bccells[0]->extrapolateQuant(EnthalpyFlow)
      +0.5*out1*bccells[1]->extrapolateQuant(EnthalpyFlow)
      -out0*bccells[0]->mHbc()();
      error.subvec(Ns*nr,Ns*(nr+1)-1)=errormHint; nr++;
    }
    {
      vd errorQ=out0*bccells[0]->extrapolateQuant(HeatFlow)
        +out1*bccells[1]->extrapolateQuant(HeatFlow);

      error.subvec(Ns*nr,Ns*(nr+1)-1)=errorQ; nr++;
    }
    {
      // d Sfgem=0.5*(bccells[0]->Sfbc()+bccells[1]->Sfbc());
      vd errorp=(bccells[0]->extrapolateQuant(Pressure)
                 -bccells[1]->extrapolateQuant(Pressure));

      error.subvec(Ns*nr,Ns*(nr+1)-1)=errorp; nr++;
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
    {                           // Mass flow continuity
      JacRow mjac(eqnr,2);
      eqnr+=Ns;
      mjac+=JacCol(bccells[0]->mbc(),out0*eye);
      mjac+=JacCol(bccells[1]->mbc(),out1*eye);
      jac+=mjac;
    }
    {
      JacRow mjacint(eqnr,5);
      eqnr+=Ns;

      // // Interpolation: mass flow is average of extrapolated mf
      mjacint+=(bccells[0]->dExtrapolateQuant(MassFlow)*=0.5*out0);
      mjacint+=(bccells[1]->dExtrapolateQuant(MassFlow)*=0.5*out1);
      mjacint+=JacCol(bccells[0]->mbc(),-out0*eye);

      jac+=mjacint;
    }
    {
      JacRow mHjac(eqnr,2);
      eqnr+=Ns;
      mHjac+=JacCol(bccells[0]->mHbc(),out0*eye);
      mHjac+=JacCol(bccells[1]->mHbc(),out1*eye);
      
      jac+=mHjac;
    }
    {
      JacRow Tjac(eqnr,2);
      eqnr+=Ns;
      Tjac+=JacCol(bccells[0]->Tbc(),eye);
      Tjac+=JacCol(bccells[1]->Tbc(),-eye);
      jac+=Tjac;

    }    
    {
      JacRow Qjac(eqnr,5);
      eqnr+=Ns;

      Qjac+=(bccells[0]->dExtrapolateQuant(HeatFlow)*=out0);
      Qjac+=(bccells[1]->dExtrapolateQuant(HeatFlow)*=out1);
      jac+=Qjac;
    }
    {
      JacRow pjac(eqnr,2);
      // eqnr+=Ns; // Not needed no row below
      pjac+=bccells[0]->dExtrapolateQuant(Pressure);
      pjac+=(bccells[1]->dExtrapolateQuant(Pressure)*=-1);
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



