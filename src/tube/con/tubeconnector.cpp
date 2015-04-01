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

namespace tube {
  using tasystem::TaSystem;

  // Number of equations corresponding to this connection:
  const int Neq=8;
  inline const RightCell& rtv(const Tube& t){
    return static_cast<const RightCell&>(t.rightCell());
  }
  inline const LeftCell& ltv(const Tube& t){
    return static_cast<const LeftCell&>(t.leftCell());
  }


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
  SimpleTubeConnector::SimpleTubeConnector(const SimpleTubeConnector& other):
    Connector(other),
    segnrs(other.segnrs),
    pos(other.pos)
  {
    TRACE(15,"SimpleTubeConnector::SimpleTubeConnector(copy)"); 
  }
  void SimpleTubeConnector::init(const TaSystem& sys){
    TRACE(15,"SimpleTubeConnector::init()");
    Connector::init(sys);
    tubes[0]=&sys.getTube(segnrs[0]);
    tubes[1]=&sys.getTube(segnrs[1]);
    setInit(true);
  }
  vd SimpleTubeConnector::error() const{
    TRACE(10,"SimpleTubeConnector::error()");

    vd error(Neq*gc->Ns());
    vd errorMf(gc->Ns(),fillwith::zeros);
    vd errorMom(gc->Ns(),fillwith::zeros);
    if(pos[0]==Pos::right){
      errorMf+=tubes[0]->rightCell().extrapolateQuant(massFlow);
    }
    else{
      errorMf-=tubes[0]->leftCell().extrapolateQuant(massFlow);
    }
    if(pos[1]==Pos::right){
      errorMf+=tubes[1]->rightCell().continuity().massFlow();
    }
    if(pos[1]==Pos::left){
      errorMf-=tubes[1]->leftCell().continuity().massFlow();
    }
    error.subvec(0,gc->Ns()-1)=errorMf;
    return error;
  }
  void SimpleTubeConnector::setEqNrs(us firstdofnr){
    TRACE(15,"SimpleTubeConnector::setEqNrs()");

  }
  us SimpleTubeConnector::getNEqs() const {
    TRACE(15,"SimpleTubeConnector::getNEqs()");
    return Neq*gc->Ns();
  }
  void SimpleTubeConnector::show(us) const{
    TRACE(15,"SimpleTubeConnector::show()");
    checkInit();
  }
  void SimpleTubeConnector::jac(tasystem::Jacobian&) const {

  }
  void SimpleTubeConnector::updateNf(){
    TRACE(15,"SimpleTubeConnector::updateNf()");

  }

} // namespace tube
//////////////////////////////////////////////////////////////////////



