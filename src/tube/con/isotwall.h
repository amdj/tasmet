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
    us firsteqnr;
    PrescribeQty massflowzero;       // Equation to make volume flow zero
                                // at wall
    PrescribeQty enthalpyflowzero;       // Equation to make volume flow zero
                                // at wall
    PrescribeQty Tbc;
    PrescribeQty Tsbc;    
    bool arbitrateMass=false;

  public:
    IsoTWall(us segnr,Pos position,const variable::var& Tbc,bool arbitrateMass=false);
  private:
    IsoTWall(const IsoTWall& o,const tasystem::TaSystem&);
  public:    
    IsoTWall(const IsoTWall& o)=delete;
    IsoTWall& operator=(const IsoTWall&) =delete;
    virtual ~IsoTWall(){}
    int arbitrateMassEq() const;
    virtual segment::Connector* copy(const tasystem::TaSystem& sys) const {
      return new IsoTWall(*this,sys);}
    #ifndef SWIG
    us getNEqs() const;
    virtual void updateNf();
    virtual void setEqNrs(us firstdofnr);    
    virtual vd error() const;
    virtual void jac(tasystem::Jacobian&) const;
    // ------------------------------
    virtual void show(us i) const;
    #endif
  };

} // namespace tube


#endif /* _ISOTWALL_H_ */




