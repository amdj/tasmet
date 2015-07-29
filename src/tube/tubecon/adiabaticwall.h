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
    // Is this connector going to arbitrate mass, yes or no?
    PrescribeQty massflowzero;
    PrescribeQty enthalpyflowzero;
    bool arbitrateMass;
  public:
    AdiabaticWall(us segnr,Pos position,bool arbitrateMass=false);
  private:
    AdiabaticWall(const AdiabaticWall& o,const tasystem::TaSystem& sys);
  public:
    AdiabaticWall& operator=(const AdiabaticWall&)=delete;
    AdiabaticWall(const AdiabaticWall& o)=delete;

    int arbitrateMassEq() const;
    virtual ~AdiabaticWall(){}
    virtual segment::Connector* copy(const tasystem::TaSystem&) const;
    virtual vd error() const;

    #ifndef SWIG
    us getNEqs() const;
    virtual void updateNf();
    virtual void setEqNrs(us firstdofnr);    
    virtual void jac(tasystem::Jacobian&) const;
    // ------------------------------
    virtual void show(us i) const;
    #endif
  };

} // namespace tube


#endif /* _ADIABATICWALL_H_ */




