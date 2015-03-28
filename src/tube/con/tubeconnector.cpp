// tubeconnector.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
// Connect two tube-type segments by using conservation of mass,
// momentum and energy.
//////////////////////////////////////////////////////////////////////

#include "tubeconnector.h"
#include "tasystem.h"
#include "tube.h"
#include "tubebcvertex.h"
#include "constants.h"

namespace tube {
  using tasystem::TaSystem;

  SimpleTubeConnector::SimpleTubeConnector(us seg1,Pos pos1,\
                                           us seg2,Pos pos2)
    
  {
    TRACE(15,"SimpleTubeConnector::SimpleTubeConnector()");
  
    if(max(seg1,seg2)>constants::maxsegs)
      throw MyError("Too high segment number given");
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

    vd error(gc->Ns());
    vd errorMf(gc->Ns(),fillwith::zeros);
    if(pos[0]==Pos::right){
      errorMf+=tubes[0]->rightVertex().extrapolateQuant(massFlow);
    }
    else{
      // errorMf+=tubes[0]->leftVertex()->massFlow();
    }
    return error;
  }
  void SimpleTubeConnector::setEqNrs(us firstdofnr){
    TRACE(15,"SimpleTubeConnector::setEqNrs()");

  }
  us SimpleTubeConnector::getNEqs() const {
    TRACE(15,"SimpleTubeConnector::getNEqs()");

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



