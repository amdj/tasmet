#pragma once
#ifndef _ADIABATICWALL_H_
#define _ADIABATICWALL_H__

#include "tubebc.h"
#include "prescribeqty.h"

namespace tasystem{
  class TaSystem;
}
namespace utils{
  template<typename SegType,typename Sys>
  SegType* copySeg(const SegType& t,const Sys& sys);
}
namespace tube{

 // Adiabatic wall boundary
  class AdiabaticWall:public TubeBc {
    bool isentropic=false;
    us firsteqnr;
    PrescribeQty massflowzero;
    PrescribeQty enthalpyflowzero;
    AdiabaticWall(const AdiabaticWall& o)=delete;
  public:
    AdiabaticWall& operator=(const AdiabaticWall&)=delete;
    AdiabaticWall(us segnr,Pos position): TubeBc(segnr,position){}
    AdiabaticWall(const AdiabaticWall& o,const tasystem::TaSystem& sys);
    virtual ~AdiabaticWall(){}
    virtual segment::Connector* copy(const tasystem::TaSystem&) const;

    #ifndef SWIG
    us getNEqs() const;
    virtual vd error() const;
    virtual void updateNf();
    virtual void setEqNrs(us firstdofnr);    
    virtual void jac(tasystem::Jacobian&) const;
    // ------------------------------
    virtual void show(us i) const;
    #endif
  };

} // namespace tube


#endif /* _ADIABATICWALL_H_ */




