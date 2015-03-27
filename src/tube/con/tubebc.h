#pragma once
#ifndef _TUBEBC_H_
#define _TUBEBC_H_
#include "connector.h"
#include "pos.h"

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
    #ifndef SWIG
    us getNEqs() const;         // 4 times Ns
    virtual ~TubeBc(){}
    virtual bool init(const tasystem::TaSystem&);
    #endif
  };

} // namespace tube



#endif /* _TUBEBC_H_ */

