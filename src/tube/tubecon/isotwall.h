#pragma once
#ifndef _ISOTWALL_H_
#define _ISOTWALL_H__

#include "tubebc.h"
#include "pressurebc.h"         // for coldtemp function
#include "prescribeqty.h"
#include "var.h"

#ifndef SWIG
namespace tasystem{
  class TaSystem;
}
#endif

namespace tube{

 // IsoT wall boundary
  class IsoTWall:public TubeBc {
    // Equation to make volume flow zero at wall
    PrescribeQty massflowzero;
    // Equation to make volume flow zero at wall
    PrescribeQty enthalpyflowzero;

    // 
    PrescribeQty Tbc;
    // Not yet implemented!
    PrescribeQty Tsbc;
    
    // 
    bool arbitrateMass=false;
  private:
    IsoTWall(const IsoTWall& o,const tasystem::TaSystem&);
  public:
    // segnr: segment number to apply this b.c. to. Position: left or
    // right side.
    // Tbc: temperature boundary condition value.
    // arbitrateMass: whether this b.c. returns an equation number
    // such that a TaSystem can overwrite this equation with global
    // mass conservation.
    // provide
    IsoTWall(us segnr,Pos position,const tasystem::var& Tbc
             ,bool arbitrateMass=false);
    IsoTWall(const IsoTWall& o)=delete;
    virtual segment::Connector* copy(const tasystem::TaSystem& sys) const {
      return new IsoTWall(*this,sys);}
    IsoTWall& operator=(const IsoTWall&) =delete;
    virtual ~IsoTWall(){}
    int arbitrateMassEq() const;

    virtual vd error() const;
    #ifndef SWIG
    us getNEqs() const;
    virtual void updateNf();
    virtual void setEqNrs(us firstdofnr);    
    virtual void jac(tasystem::Jacobian&) const;
    virtual void show(us i) const;
    #endif
  };

} // namespace tube


#endif /* _ISOTWALL_H_ */




