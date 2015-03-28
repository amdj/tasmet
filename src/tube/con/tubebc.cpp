#include "tubebc.h"
#include "tasystem.h"
#include "tube.h"

namespace tube{
  using tasystem::TaSystem;

  void TubeBc::init(const TaSystem& sys){
    TRACE(15,"TubeBc::init()");
    this->sys=&sys;
    Connector::init(sys);
    try{
      t=&sys.getTube(segnr);
    }
    catch(std::exception& b){
      WARN("Seg nr " << segnr << " is not a Tube! Initialization failed.");
      WARN("In: " << getType() << " with name "<< getName() <<", and number:" << getNumber());
      throw;
    }
  } // init
  us TubeBc::getNEqs() const {
      TRACE(10,"TubeBc::getNEqs()");
      return 4*gc->Ns();
    }

} // namespace tube
