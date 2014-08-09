#include "tubebcvertex.h"

namespace tube{

  TubeBcVertex::TubeBcVertex(us segnr){
    segnumber=segnr;  
  }
  TubeBcVertex::TubeBcVertex(const TubeBcVertex& o):
    TubeBcVertex(o.segNumber())
  {
    TRACE(8,"TubeBcVertex copy cc.");
  }
  TubeBcVertex& TubeBcVertex::operator=(const TubeBcVertex& o)
    {
      TRACE(8,"TubeBcVertex::operator=()");
      setSegNumber(o.segNumber());
      return *this;      
    }
    // virtual Vertex* copy(const SegBase&)=0; // Copy the boundary condition vertex
  
} // namespace tube


