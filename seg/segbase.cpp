#include "segbase.h"

namespace segment{

  SegBase::SegBase(const Geom& geom):geom_(geom.copy())
  {
    TRACE(10,"SegBase::SegBase(geom)");
  }
  SegBase::SegBase(const SegBase& o): SegBase(o.geom()){}
  void SegBase::init(const Globalconf& gc1){this->gc=&gc1;}  
  SegBase::~SegBase(){
    delete geom_;
  }


} // namespace segment
