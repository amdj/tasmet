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
#include "leftcell.h"
#include "rightcell.h"

#define fDFT (gc->fDFT)
#define iDFT (gc->iDFT)
#define DDTfd (gc->DDTfd)

#define Ns (gc->Ns())
#define eye (eye(Ns,Ns))

namespace tube {
  using tasystem::TaSystem;

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
    tubes[0]=&sys.getTube(segnrs[0]);
    tubes[1]=&sys.getTube(segnrs[1]);
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
    vd errorMf(Ns,fillwith::zeros);
    vd errorMfint(Ns,fillwith::zeros);    
    vd errorMom(Ns,fillwith::zeros);

    errorMfint+=out0*tubes[0]->bcCell(pos[0]).extrapolateQuant(massFlow);
    errorMfint+=out1*tubes[1]->bcCell(pos[1]).extrapolateQuant(massFlow);

    if(pos[1]==Pos::right){
      // errorMf+=tubes[1]->rightCell().massFlow();
    }
    if(pos[1]==Pos::left){
      // errorMf-=tubes[1]->leftCell().continuity().massFlow();
    }
    error.subvec(0,Ns-1)=errorMfint;
    return error;
  }
  void SimpleTubeConnector::setEqNrs(us firstdofnr){
    TRACE(15,"SimpleTubeConnector::setEqNrs()");

  }
  us SimpleTubeConnector::getNEqs() const {
    TRACE(15,"SimpleTubeConnector::getNEqs()");
    return Neq*Ns;
  }
  void SimpleTubeConnector::show(us) const{
    TRACE(15,"SimpleTubeConnector::show()");
    checkInit();
  }
  void SimpleTubeConnector::jac(tasystem::Jacobian&) const {
    TRACE(15,"SimpleTubeconnector::jac()");
    WARN("Does nothing");
  }
  void SimpleTubeConnector::updateNf(){
    TRACE(15,"SimpleTubeConnector::updateNf()");
  }

} // namespace tube
//////////////////////////////////////////////////////////////////////



