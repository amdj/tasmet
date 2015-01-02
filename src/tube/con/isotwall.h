#pragma once
#ifndef _ISOTWALL_H_
#define _ISOTWALL_H__

#include "tubebc.h"
#include "pressurebc.h"         // for coldtemp function
#include "prescribeqty.h"
#include "var.h"

namespace tasystem{
  class TaSystem;
}

namespace tube{

 // IsoT wall boundary
  class IsoTWall:public TubeBc {
    us firsteqnr;
    vd zero;
    PrescribeQty Uiszero;       // Equation to make volume flow zero
                                // at wall
    PrescribeQty Tbc;
    PrescribeQty Tsbc;    
    IsoTWall& operator=(const IsoTWall&);
  public:
    IsoTWall(us segnr,pos position,const variable::var& Tbc,const variable::var& Tsbc);
    IsoTWall(us segnr,pos position,const variable::var& Tbc):
      IsoTWall(segnr,position,Tbc,coldtemp(Tbc)){}
    IsoTWall(const IsoTWall& o);
    virtual ~IsoTWall(){}

    virtual segment::Connector* copy() const {return new IsoTWall(*this);}
    virtual string getType() const {return string("IsoTWall");}

    virtual bool init(const tasystem::TaSystem&);
    virtual void updateNf();
    virtual void setEqNrs(us firstdofnr);    
    virtual vd error() const;
    virtual void jac(tasystem::Jacobian&) const;
    // ------------------------------
    virtual void show(us i) const;

  };

} // namespace tube


#endif /* _ISOTWALL_H_ */




