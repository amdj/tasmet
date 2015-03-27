#include "segconbase.h"
#include "tasystem.h"

namespace segment{
  using tasystem::Globalconf;   
  using tasystem::TaSystem;   


  us SegConBase::globnr_=0;
  SegConBase::SegConBase(){
    TRACE(13,"SegConBase::SegConBase()");
    globnr_++;
    string s="Nameless";
    s+=std::to_string(globnr_);
    name_=s;
  }
  SegConBase::SegConBase(const SegConBase& o):
    name_(o.name_)
  {/* Copy constructor */  }
  bool SegConBase::init(const TaSystem& sys){
    TRACE(13,"Seg::init()");
    this->gc=&sys.gc;
    return true;
  }

} // namespace segment
