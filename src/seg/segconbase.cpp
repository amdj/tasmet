#include "segconbase.h"
#include "tasystem.h"
#include <typeinfo>

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
  SegConBase::SegConBase(const SegConBase& o,const TaSystem& sys):
    name_(o.name_)
  {/* Copy constructor */
    this->gc=&sys.gc();
  }
  void SegConBase::setNumber(us number) {
    TRACE(15,"setNumber called for type " << typeid(*this).name());
      this->number=number;} 

} // namespace segment
