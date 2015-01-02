#pragma once
#ifndef _ADIABATICWALL_H_
#define _ADIABATICWALL_H__

#include "tubebc.h"
#include "prescribeqty.h"

namespace tasystem{
  class TaSystem;
}

namespace tube{

 // Adiabatic wall boundary
  class AdiabaticWall:public TubeBc {
    vd zero;
    PrescribeQty Uiszero;       // Equation to make volume flow zero
                                // at wall
    PrescribeddxQty drhodxiszero;
    us firsteqnr;
    AdiabaticWall& operator=(const AdiabaticWall&);
  public:
    AdiabaticWall(us segnr,pos position): TubeBc(segnr,position){}
    AdiabaticWall(const AdiabaticWall& o): TubeBc(o) {}
    virtual ~AdiabaticWall(){}

    virtual segment::Connector* copy() const {return new AdiabaticWall(*this);}
    virtual string getType() const {return string("AdiabaticWall");}

    virtual bool init(const tasystem::TaSystem&);
    virtual void updateNf();
    virtual void setEqNrs(us firstdofnr);    
    virtual vd error() const;
    virtual void jac(tasystem::Jacobian&) const;
    // ------------------------------
    virtual void show(us i) const;

  };

} // namespace tube


#endif /* _ADIABATICWALL_H_ */




