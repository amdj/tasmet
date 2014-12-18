#include "segconbase.h"
#include "tasystem.h"

namespace segment{
  using tasystem::Globalconf;   
  using tasystem::TaSystem;   


  us SegConBase::globnr_=0;
  SegConBase::SegConBase(){
    TRACE(13,"SegConBase::SegConBase()");
    globnr_++;
    std::stringstream s;
    s << "Nameless" << globnr_;
    name_=s.str();
  }
  SegConBase::SegConBase(const SegConBase& o):
    name_(o.name_)
  {/* Copy constructor */  }
  void SegConBase::init(const TaSystem& sys){
    TRACE(13,"Seg::init()");
    this->gc=&sys.gc;
    init_=true;
  }  

} // namespace segment
