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
  protected:
    Connector():SegConBase(){}
    #ifndef SWIG
    Connector& operator=(const Connector&)=delete;
    Connector(const Connector&)=delete;
    #endif // ifndef SWIG
    Connector(const Connector& o,const tasystem::TaSystem& sys):SegConBase(o,sys){}
  public:
    virtual ~Connector(){}
    virtual Connector* copy(const tasystem::TaSystem&) const=0;
    #ifndef SWIG
    // 100% Forwarding

    #endif
  };

} // namespace segment

#endif /* _CONNECTOR_H_ */
