#pragma once
#ifndef _CONNECTOR_H_
#define _CONNECTOR_H_
#include "seg.h"

namespace variable{class var;}

namespace segment{

  class Connector:public Seg{
  public:
    const variable::var& rho() const;
    const variable::var& U()   const;
    const variable::var& p()   const;
    const variable::var& T()   const;
    const variable::var& Ts()  const;
  };

} // namespace segment

#endif /* _CONNECTOR_H_ */
