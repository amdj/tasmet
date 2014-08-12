// File vertex.h
#pragma once
#ifndef _VERTEX_H_
#define _VERTEX_H_

#include "globalconf.h"
#include "var.h"


namespace segment{
  SPOILNAMESPACE
  class Seg;
  using variable::var;
  using tasystem::Globalconf;
  
  class Vertex;
  
  typedef std::vector<const Vertex*>  VertexVec;
  class Vertex{
  public:
    const tasystem::Globalconf* gc=NULL;
    us i;

    virtual ~Vertex(){}
    // const VertexVec& getLeft() const {return vleft;}
    // const VertexVec& getRight() const {return vright;}
    virtual vd error() const=0;		       // Compute error for this gridpoint
    virtual dmat jac() const=0;		       // Fill complete Jacobian for this node
    virtual void setRes(vd res)=0;			  // Set result vector to res
    virtual vd getRes() const=0;			  // Extract current result vector
    virtual void show() const=0;
    Vertex& operator=(const Vertex& v2){WARN("Operator=() not allowed. Aborting."); abort();} 
    void init(us i,const Globalconf& gc);
  };
} // namespace segment


#endif /* _VERTEX_H_ */
