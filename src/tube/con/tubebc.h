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
  #ifdef SWIG
  %feature("abstract") TubeBc;
  #endif
  class TubeBc:public segment::Connector {
  protected:
    us segnr;
    Pos pos;
    const Tube* t=nullptr;
    const tasystem::TaSystem* sys=nullptr;
  public:
    TubeBc(us segnr,Pos position):segnr(segnr),pos(position){}
    TubeBc(const TubeBc& other):
      Connector(other),
      segnr(other.segnr),
      pos(other.pos)
    {}
    virtual ~TubeBc(){}
    #ifndef SWIG
    us getNEqs() const;         // 4 times Ns
    virtual void init(const tasystem::TaSystem&);
    #endif
  };

} // namespace tube



#endif /* _TUBEBC_H_ */

