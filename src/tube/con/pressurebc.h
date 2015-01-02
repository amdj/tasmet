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

  variable::var coldtemp(const variable::var&);

  class PressureBc:public TubeBc {
    us firsteqnr;
    PrescribeQty prescribep;			// Pressure boundary condition
    PrescribeQty prescribeT;			// Temperature boundary condition
    PrescribeQty prescribeTs;			// Solid temperature boundary condition
    
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
    virtual bool init(const tasystem::TaSystem&);
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

  };

} // namespace tube

#endif /* _PRESSUREBC_H_ */


