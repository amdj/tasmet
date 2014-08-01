// File vertex.h
#pragma once
#ifndef _VERTEX_H_
#define _VERTEX_H_

#include "segbase.h"
#include "var.h"
#include "localgeom.h"


namespace segment{
  SPOILNAMESPACE
  class Seg;
  using segment::Geom;
  using variable::var;

  
  class Vertex{
  public:
    LocalGeom lg;
    const Globalconf* gc=NULL;
    Vertex(){}
    virtual vd error() const=0;		       // Compute error for this gridpoint
    virtual dmat jac() const=0;		       // Fill complete Jacobian for this node
    virtual void setRes(vd res)=0;			  // Set result vector to res
    virtual vd getRes() const=0;			  // Extract current result vector
    virtual void show() const=0;
    virtual ~Vertex() {}
    Vertex(const Vertex&); 
    Vertex& operator=(const Vertex& v2){WARN("Assigning operators not allowed. Aborting."); abort();} 
    void init(us i,const SegBase& thisseg);

  };
} // namespace segment


#endif /* _VERTEX_H_ */
