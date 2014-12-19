#pragma once
#ifndef _PRESSUREBC_H_
#define _PRESSUREBC_H_
#include "var.h"
#include "tubebc.h"
#include "prescribeqty.h"

namespace tasystem{
  class TaSystem;
}

namespace tube{
  class PressureBc:public TubeBc {
    us firsteqnr;
    variable::var pbc;			// Pressure boundary condition
    variable::var Tbc;			// Temperature boundary condition
    variable::var Tsbc;			// Solid temperature boundary condition
    PrescribeQty prescribep;
    PrescribeQty prescribeT;
    PrescribeQty prescribeTs;
    
    PressureBc& operator=(const PressureBc&);
  public:
    // Set all variables
    PressureBc(const variable::var& p,const variable::var& T,const variable::var& Ts,us segnr,pos position);
    // Assume solid temperature constant at gc.T0;
    PressureBc(const variable::var& p,const variable::var& T,us segnr,pos position); 
    // Assume above and adiabatic compresion/expansion
    PressureBc(const variable::var& p,us segnr,pos position);
    PressureBc(const PressureBc& other);
    virtual ~PressureBc(){}
    virtual void init(const tasystem::TaSystem&);
    virtual segment::Connector* copy() const { return new PressureBc(*this);}
    virtual string getType() const {return string("PressureBc");}

    virtual void updateNf();
    virtual void setEqNrs(us firstdofnr);    
    virtual vd error() const;
    virtual void jac(tasystem::Jacobian&) const;
    // ------------------------------
    virtual void show(us i) const;
  private:
    static variable::var adiabatictemp(const variable::var& pres); // Return adiabatic compression
    // amplitude values
    static variable::var coldtemp(const variable::var& pres); // Returns gc.T0 amplitude data

  };

} // namespace tube

#endif /* _PRESSUREBC_H_ */


