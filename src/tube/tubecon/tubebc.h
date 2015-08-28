#pragma once
#ifndef _TUBEBC_H_
#define _TUBEBC_H_
#include "connector.h"
#include "constants.h"

namespace tube{

  #ifndef SWIG
  class Tube;
  #endif

  class TubeBc:public segment::Connector {
  protected:
    std::string segid;
    Pos pos;
    us firsteqnr;
  protected:
    const Tube* t=nullptr;
    const tasystem::TaSystem* sys=nullptr;
    TubeBc(const TubeBc& other,const tasystem::TaSystem& sys);
    TubeBc(const string& segid,Pos position):segid(segid),pos(position){}
    void setEqNrs(us firsteqnr1) { firsteqnr=firsteqnr1;}
  public:
    virtual ~TubeBc(){}
    #ifndef SWIG
    us getNEqs() const=0;
    #endif
  };

} // namespace tube



#endif /* _TUBEBC_H_ */

