#include "segbase.h"

namespace segment{
  static us totalnumber=0;

  SegBase::SegBase(Geom geom):geomptr(new Geom(geom)),geom(*geomptr)
  {
    number=totalnumber;
    nleft=0;
    nright=0;
    totalnumber++;
 
  }
  SegBase::~SegBase(){
    TRACE(-5,"~SegBase()");
    delete geomptr;
  }
  void SegBase::newgeom(const Geom& g){
    delete geomptr;
    geomptr=new Geom(g);
  }
  const string& SegBase::gettype() const {return type;}
  bool SegBase::operator==(const SegBase& other) const {return (this->number==other.number);}


} // namespace segment
