#include "tubebc.h"
#include "tasystem.h"
#include "tube.h"

namespace tube{
  using tasystem::TaSystem;

  bool TubeBc::init(const TaSystem& sys){
    TRACE(15,"TubeBc::init()");
    this->sys=&sys;
    if(!Connector::init(sys))
      return false;
    try{
      t=dynamic_cast<const Tube*>(sys.getSeg(segnr));
    }
    catch(std::bad_cast){
      WARN("Seg nr " << segnr << " is not a Tube! Initialization failed.");
      return false;
    }
    return true;
  } // init
  us TubeBc::getNEqs() const {
      TRACE(10,"TubeBc::getNEqs()");
      return 4*gc->Ns();
    }

} // namespace tube
