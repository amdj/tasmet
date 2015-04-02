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
    bool isentropic=false;
    us firsteqnr;
    PrescribeQty massflowzero;
  public:
    AdiabaticWall& operator=(const AdiabaticWall&)=delete;
    AdiabaticWall(us segnr,Pos position): TubeBc(segnr,position){}
    AdiabaticWall(const AdiabaticWall& o): TubeBc(o) {}
    virtual ~AdiabaticWall(){}
    virtual segment::Connector* copy() const {return new AdiabaticWall(*this);}
    virtual string getType() const {return string("AdiabaticWall");}

    #ifndef SWIG
    us getNEqs() const;
    virtual vd error() const;
    virtual void init(const tasystem::TaSystem&);
    virtual void updateNf();
    virtual void setEqNrs(us firstdofnr);    
    virtual void jac(tasystem::Jacobian&) const;
    // ------------------------------
    virtual void show(us i) const;
    #endif
  };

} // namespace tube


#endif /* _ADIABATICWALL_H_ */




