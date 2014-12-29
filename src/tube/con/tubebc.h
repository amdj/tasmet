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
  protected:
    us segnr;
    pos position;
    const Tube* t=nullptr;
    const tasystem::TaSystem* sys=nullptr;
  public:
    TubeBc(us segnr,pos position):segnr(segnr),position(position){}
    TubeBc(const TubeBc& other):
      Connector(other),
      segnr(other.segnr),
      position(other.position)
    {}
    us getNEqs() const;         // 4 times Ns
    virtual ~TubeBc(){}
    virtual bool init(const tasystem::TaSystem&);
  };

} // namespace tube



#endif /* _TUBEBC_H_ */

