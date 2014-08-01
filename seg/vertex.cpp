#include "vertex.h"
#include "var.h"


namespace segment{

  Vertex::Vertex(const Vertex& o):Vertex(){
    TRACE(8,"Vertex copy constructor. Copying only LocalGeom");
    lg=o.lg;
  }
  void Vertex::init(us i,const SegBase& thisseg){
    TRACE(8,"Vertex::Init()");
    gc=thisseg.gc;
    lg=thisseg.geom.localGeom(i);
  }


} // namespace segment

