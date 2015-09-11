#pragma once
#ifndef _DUCTBC_H_
#define _DUCTBC_H_
#include "connector.h"
#include "constants.h"

namespace duct{

  #ifndef SWIG
  class Duct;
  #endif

  class DuctBc:public segment::Connector {
  protected:
    std::string segid;
    Pos pos;
    us firsteqnr;
  protected:
    const Duct* t=nullptr;
    const tasystem::TaSystem* sys=nullptr;
    DuctBc(const DuctBc& other,const tasystem::TaSystem& sys);
    DuctBc(const string& segid,Pos position):segid(segid),pos(position){}
    void setEqNrs(us firsteqnr1) { firsteqnr=firsteqnr1;}
  public:
    virtual ~DuctBc(){}
    #ifndef SWIG
    us getNEqs() const=0;
    #endif
  };

} // namespace duct



#endif /* _DUCTBC_H_ */

