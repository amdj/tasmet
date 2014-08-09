#include "segbase.h"

namespace segment{

  SegBase::SegBase(const Geom& geom):geom(geom)
  {
    TRACE(10,"SegBase::SegBase(geom)");
  }
  SegBase::SegBase(const SegBase& o): SegBase(o.geom){}
  SegBase& SegBase::operator=(const SegBase& o){geom=o.geom; return *this;}
  void SegBase::init(const Globalconf& gc1){this->gc=&gc1;}  

  void SegBase::setLeft(const SegBase& Left){
    TRACE(13,"SegBase::SetLeft()");
    left.push_back(&Left);
  }
  void SegBase::setRight(const SegBase& Right){
    TRACE(13,"SegBase::SetRight()");
    right.push_back(&Right);
  }
  

  bool SegBase::operator==(const SegBase& other) const {return (this->number==other.number);}


} // namespace segment
