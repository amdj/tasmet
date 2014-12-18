#pragma once
#ifndef _CONNECTOR_H_
#define _CONNECTOR_H_
#include "vtypes.h"
#include "segconbase.h"

namespace tasystem{
  class TaSystem;
}
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
    // 100% Forwarding
    virtual void init(const tasystem::TaSystem& sys){SegConBase::init(sys);}
  };

} // namespace segment

#endif /* _CONNECTOR_H_ */
