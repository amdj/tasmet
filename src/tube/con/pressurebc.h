#pragma once
#ifndef _PRESSUREBC_H_
#define _PRESSUREBC_H_
#include "var.h"
#include "connector.h"
#include "pos.h"

namespace tasystem{
  class TaSystem;
}

namespace tube{
  class PressureBc:public segment::Connector
  {
    us segnr;
    pos position;
    variable::var pLbc;			// Pressure boundary condition
    variable::var TLbc;			// Temperature boundary conditions
    PressureBc& operator=(const PressureBc&);
  public:
    PressureBc(const variable::var& p,const variable::var& T,us segnr,pos position);
    PressureBc(const variable::var& p,us segnr,pos position);
    PressureBc(const PressureBc& other);
    virtual ~PressureBc(){}
    virtual void init(const tasystem::TaSystem&);
    virtual Connector* copy() const { return new PressureBc(*this);}
    virtual string getType() const {return string("PressureBc");}

    virtual void updateNf();
    // ------------------------------
    virtual void setEqNrs(us firstdofnr);    
    virtual void show() const;
    virtual us getNEqs() const;    
    virtual vd error() const;
    virtual void show(us) const;
    virtual void jac(tasystem::Jacobian&) const;
  };

} // namespace tube

#endif /* _PRESSUREBC_H_ */


