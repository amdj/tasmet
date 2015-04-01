#include "segconbase.h"
#include "tasystem.h"

namespace segment{
  using tasystem::Globalconf;   
  using tasystem::TaSystem;   


  us SegConBase::globnr_=0;
  SegConBase::SegConBase(){
    TRACE(13,"SegConBase::SegConBase()");
    globnr_++;
    string s="Nameless ";
    s+=std::to_string(globnr_);
    name_=s;
  }
  SegConBase::SegConBase(const SegConBase& o):
    name_(o.name_)
  {/* Copy constructor */
  }
  void SegConBase::init(const TaSystem& sys){
    TRACE(13,"SegConBase::init()");
    this->gc=&sys.gc();
  }
  void SegConBase::setNumber(us number) {
      TRACE(15,"setNumber called for type" << getType());
      this->number=number;} 

} // namespace segment
