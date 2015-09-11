#include "ductbc.h"
#include "tasystem.h"
#include "duct.h"
#include "exception.h"

namespace duct{
  using tasystem::TaSystem;


  DuctBc::DuctBc(const DuctBc& other,const TaSystem& sys):
      Connector(other,sys),
      segid(other.segid),
      pos(other.pos)
  {
    this->sys=&sys;
    try{
      t=&sys.getDuct(segid);
    }
    catch(const std::exception& b){
      WARN("Seg id" << segid << " is not a Duct! Initialization of DuctBc failed.");
      WARN("In: " << typeid(*this).name() << " with name "<< getName() <<", and number:" << getNumber());
      throw MyError(b.what());
    }
  } // init

} // namespace duct
