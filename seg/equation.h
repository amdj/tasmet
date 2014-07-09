#pragma once
#include <vtypes.h>
#include "globalconf.h"


namespace segment{
  SPOILNAMESPACE
  class Vertex;

  class Equation
  {
  public:
    const tasystem::Globalconf* gc=NULL;
    Vertex& vertex;
    
    Equation(Vertex& v):vertex(v) {}
    ~Equation(){}

    virtual vd Error()=0;
    virtual dmat Jac()=0;		// Returns the local Jacobian of this equation
    virtual void Init(const tasystem::Globalconf& g){gc=&g;}
  };

  
}




