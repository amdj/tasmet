#include "segbase.h"

namespace segment{

  SegBase::SegBase(Geom geom):geom(geom)
  {
  }
  SegBase::~SegBase(){
    TRACE(-5,"~SegBase()");
    }
  
  void SegBase::newgeom(const Geom& g){
    // delete geomptr;
    // geomptr=new Geom(g);
  }
  void SegBase::setLeft(const SegBase& Left){
    TRACE(13,"SegBase::SetLeft()");
    left.push_back(&Left);
  }
  void SegBase::setRight(const SegBase& Right){
    TRACE(13,"SegBase::SetRight()");
    right.push_back(&Right);
  }
  
  const string& SegBase::gettype() const {return type;}
  bool SegBase::operator==(const SegBase& other) const {return (this->number==other.number);}


} // namespace segment
