#pragma once
#ifndef _TUBEBC_H_
#define _TUBEBC_H_
#include "connector.h"
#include "pos.h"

namespace tasystem{
  class TaSystem;
}

namespace tube{

  class Tube;
  
  class TubeBc:public segment::Connector {
    us segnr,firsteqnr;
    pos position;
    const Tube* t=nullptr;
    const tasystem::TaSystem* sys=nullptr;
  public:
    TubeBc(us segnr,pos position):segnr(segnr),position(position){}
    TubeBc(const TubeBc& other):TubeBc(other.segnr,other.position){}
    virtual ~TubeBc(){}
    virtual void init(const tasystem::TaSystem&);
  };

} // namespace tube



#endif /* _TUBEBC_H_ */

