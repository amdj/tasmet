#include "vertex.h"
#include "var.h"


namespace segment{

  void Vertex::init(us i,const Globalconf& gc1){
    TRACE(8,"Vertex::Init()");
    this->i=i;
    TRACE(100,"Address of gc I got:"<<&gc1);
    this->gc=&gc1;
  }
  void Vertex::setLeft(const Vertex& v) { left.push_back(&v);}
  void Vertex::setRight(const Vertex& v) {right.push_back(&v);}

} // namespace segment

