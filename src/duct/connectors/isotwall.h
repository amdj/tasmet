#pragma once
#ifndef _ISOTWALL_H_
#define _ISOTWALL_H__

#include "ductbc.h"
#include "pressurebc.h"         // for coldtemp function
#include "prescribeqty.h"
#include "var.h"

#ifndef SWIG
namespace tasystem{
  class TaSystem;
}
#endif

namespace duct{

 // IsoT wall boundary
  class IsoTWall:public DuctBc {
    // Equation to make volume flow zero at wall
    PrescribeQty massflowzero;
    // Equation to make volume flow zero at wall
    PrescribeQty enthalpyflowzero;

    PrescribeQty Tbc;
    // Not yet implemented!
    PrescribeQty Tsbc;

    bool arbitrateMass=false;
    #ifndef SWIG
    IsoTWall(const IsoTWall& o)=delete;
    IsoTWall& operator=(const IsoTWall&) =delete;
    #endif // ifndef SWIG
    IsoTWall(const IsoTWall& o,const tasystem::TaSystem&);
  public:
    // segid: segment number to apply this b.c. to. Position: left or
    // right side.
    // Tbc: temperature boundary condition value.
    // arbitrateMass: whether this b.c. returns an equation number
    // such that a TaSystem can overwrite this equation with global
    // mass conservation.
    // provide
    IsoTWall(const string& segid,Pos position,const tasystem::var& Tbc
             ,bool arbitrateMass=false);
    virtual segment::Connector* copy(const tasystem::TaSystem& sys) const {
      return new IsoTWall(*this,sys);}

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

} // namespace duct


#endif /* _ISOTWALL_H_ */




