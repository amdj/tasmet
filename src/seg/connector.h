#pragma once
#ifndef _CONNECTOR_H_
#define _CONNECTOR_H_
#include "vtypes.h"
#include "segconbase.h"


#ifndef SWIG
namespace tasystem{
  class TaSystem;
}
#endif

namespace segment{
  // A connector contains only equations, no degrees of freedom
  
  class Connector:public SegConBase{
  public:
    Connector():SegConBase(){}
    Connector(const Connector& o):
      SegConBase(o)
    {}
    virtual ~Connector(){}
    virtual Connector* copy() const=0;
    #ifndef SWIG
    // 100% Forwarding
    virtual void init(const tasystem::TaSystem& sys){ SegConBase::init(sys);}
    #endif
  };

} // namespace segment

#endif /* _CONNECTOR_H_ */
