// pressurebc.h
//
// Author: J.A. de Jong 
//
// Description:
// Pressure boundary condition for a Tube
// 
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef PRESSUREBC_H
#define PRESSUREBC_H


#include "constants.h"
#include "var.h"
#include "tubebc.h"
#include "prescribeqty.h"

namespace tasystem{
  class TaSystem;
}


namespace tube{
  #ifndef SWIG
  tasystem::var coldtemp(const tasystem::var&);
  #endif
  #ifdef SWIG
  %catches(std::exception,...) PressureBc::PressureBc(const tasystem::var& p,const tasystem::var& T,const tasystem::var& Ts,us segnr,Pos position);
  %catches(std::exception,...) PressureBc::PressureBc(const tasystem::var& p,const tasystem::var& T,us segnr,Pos position); 
  %catches(std::exception,...) PressureBc::PressureBc(const tasystem::var& p,us segnr,Pos position);  // %feature("notabstract") PressureBc;
  #endif // SWIG

  class PressureBc:public TubeBc {
    us firsteqnr;
    tasystem::var p_prescribed;
    // PrescribeQty prescribep;			// Pressure boundary condition
    PrescribeQty prescribeT;			// Temperature boundary condition
    // PrescribeQty prescribeTs;			// Solid temperature boundary condition
  public:
    PressureBc& operator=(const PressureBc&)=delete;
    // Set all variables
    PressureBc(us segnr,Pos position,const tasystem::var& p,const tasystem::var& T,const tasystem::var& Ts);
    // Assume solid temperature constant at gc.T0;
    PressureBc(us segnr,Pos position,const tasystem::var& p,const tasystem::var& T); 
    // Assume above and adiabatic compresion/expansion
    PressureBc(us segnr,Pos position,const tasystem::var& p);
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

#endif // PRESSUREBC_H
//////////////////////////////////////////////////////////////////////

