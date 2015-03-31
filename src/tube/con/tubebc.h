#pragma once
#ifndef _TUBEBC_H_
#define _TUBEBC_H_
#include "connector.h"
#include "constants.h"

#ifndef SWIG
namespace tasystem{
  class TaSystem;
}
#endif

namespace tube{

  #ifndef SWIG
  class Tube;
  #endif

  class TubeBc:public segment::Connector {
  protected:
    us segnr;
    Pos pos;
  protected:
    const Tube* t=nullptr;
    const tasystem::TaSystem* sys=nullptr;

    TubeBc(us segnr,Pos position):segnr(segnr),pos(position){}
    TubeBc(const TubeBc& other):
      Connector(other),
      segnr(other.segnr),
      pos(other.pos)
    {}
  public:
    virtual ~TubeBc(){}
    #ifndef SWIG
    us getNEqs() const;
    virtual void init(const tasystem::TaSystem&);
    #endif
  };

} // namespace tube



#endif /* _TUBEBC_H_ */

