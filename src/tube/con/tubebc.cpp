#include "tubebc.h"

namespace tube{
  using tasystem::TaSystem;

  void TubeBc::init(const TaSystem& sys){
    TRACE(15,"TubeBc::init()");
    this->sys=&sys;
    Connector::init(sys);
  }

} // namespace tube
