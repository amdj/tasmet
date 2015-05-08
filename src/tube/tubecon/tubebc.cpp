#include "tubebc.h"
#include "tasystem.h"
#include "tube.h"
#include "exception.h"

namespace tube{
  using tasystem::TaSystem;


  TubeBc::TubeBc(const TubeBc& other,const TaSystem& sys):
      Connector(other,sys),
      segnr(other.segnr),
      pos(other.pos)
  {
    this->sys=&sys;
    try{
      t=&sys.getTube(segnr);
    }
    catch(const std::exception& b){
      WARN("Seg nr " << segnr << " is not a Tube! Initialization of TubeBc failed.");
      WARN("In: " << typeid(*this).name() << " with name "<< getName() <<", and number:" << getNumber());
      throw MyError(b.what());
    }
  } // init

} // namespace tube
