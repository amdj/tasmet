#include "vertex.h"
#include "var.h"


namespace segment{

  void Vertex::init(us i,const Globalconf& gc1){
    TRACE(8,"Vertex::Init()");
    this->i=i;
    this->gc=&gc1;
  }
  // void Vertex::setLeft(const Vertex& v) { vleft.push_back(&v);}
  // void Vertex::setRight(const Vertex& v) {vright.push_back(&v);}

} // namespace segment

