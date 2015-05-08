#pragma once
#ifndef _PRESSUREBC_H_
#define _PRESSUREBC_H_
#include "constants.h"
#include "var.h"
#include "tubebc.h"
#include "prescribeqty.h"

namespace tasystem{
  class TaSystem;
}



namespace tube{
  #ifndef SWIG
  variable::var coldtemp(const variable::var&);
  #endif
  #ifdef SWIG
  %catches(std::exception,...) PressureBc::PressureBc(const variable::var& p,const variable::var& T,const variable::var& Ts,us segnr,Pos position);
  %catches(std::exception,...) PressureBc::PressureBc(const variable::var& p,const variable::var& T,us segnr,Pos position); 
  %catches(std::exception,...) PressureBc::PressureBc(const variable::var& p,us segnr,Pos position);  // %feature("notabstract") PressureBc;
  #endif // SWIG

  class PressureBc:public TubeBc {
    us firsteqnr;
    variable::var p_prescribed;
    // PrescribeQty prescribep;			// Pressure boundary condition
    PrescribeQty prescribeT;			// Temperature boundary condition
    // PrescribeQty prescribeTs;			// Solid temperature boundary condition
  public:
    PressureBc& operator=(const PressureBc&)=delete;
    // Set all variables
    PressureBc(const variable::var& p,const variable::var& T,const variable::var& Ts,us segnr,Pos position);
    // Assume solid temperature constant at gc.T0;
    PressureBc(const variable::var& p,const variable::var& T,us segnr,Pos position); 
    // Assume above and adiabatic compresion/expansion
    PressureBc(const variable::var& p,us segnr,Pos position);
    PressureBc(const PressureBc& other)=delete;
    PressureBc(const PressureBc& other,const tasystem::TaSystem&);
    segment::Connector* copy(const tasystem::TaSystem& s) const { return new PressureBc(*this,s);}
    virtual vd error() const;
    virtual ~PressureBc(){}
    #ifndef SWIG
    us getNEqs() const {return 3*gc->Ns();}    
    virtual void updateNf();
    virtual void setEqNrs(us firstdofnr);    
    virtual void jac(tasystem::Jacobian&) const;
    // ------------------------------
    virtual void show(us i) const;
  private:
    #endif // SWIG    
  };

} // namespace tube

#endif /* _PRESSUREBC_H_ */







