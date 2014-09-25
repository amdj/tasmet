#include "vertex.h"
#include "var.h"
#include "seg.h"

namespace segment{

  void Vertex::initVertex(us i,const Seg& seg){
    TRACE(8,"Vertex::Init()");
    // Initialize the Globalconf* ptr and i (the vertex number), 
    this->i=i;
    this->gc=seg.gc;
    lg=seg.geom.localGeom(i);
  }
  // void Vertex::setLeft(const Vertex& v) { vleft.push_back(&v);}
  // void Vertex::setRight(const Vertex& v) {vright.push_back(&v);}

} // namespace segment

