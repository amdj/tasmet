#include "tubebcvertex.h"

namespace tube{

  TubeBcVertex::TubeBcVertex(){  }
  TubeBcVertex::TubeBcVertex(const TubeBcVertex& o):
    TubeBcVertex()
  {
    TRACE(8,"TubeBcVertex copy cc.");
  }
  TubeBcVertex& TubeBcVertex::operator=(const TubeBcVertex& o)
    {
      TRACE(8,"TubeBcVertex::operator=()");
      return *this;      
    }
  // virtual Vertex* copy(const SegBase&)=0; // Copy the boundary condition vertex
  
} // namespace tube


