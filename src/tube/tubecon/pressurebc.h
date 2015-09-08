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
#include "varutils.h"

namespace tasystem{
  class TaSystem;
}


namespace tube{
  #ifdef SWIG
  %catches(std::exception,...) PressureBc::PressureBc(const string& segid,Pos position,const tasystem::var& p,const tasystem::var& T,const tasystem::var& Ts);
  %catches(std::exception,...) PressureBc::PressureBc(const string& segid,Pos position,const tasystem::var& p,const tasystem::var& T); 
  %catches(std::exception,...) PressureBc::PressureBc(const string& segid,Pos position,const tasystem::var& p);  // %feature("notabstract") PressureBc;
  #endif // SWIG

  class PressureBc:public TubeBc {
    tasystem::var p_prescribed;
    // PrescribeQty prescribep;			// Pressure boundary condition
    PrescribeQty prescribeT;			// Temperature boundary condition
    PrescribeQty prescribeTs;			// Solid temperature boundary condition
    PressureBc(const PressureBc& other,const tasystem::TaSystem&);
  public:
    PressureBc& operator=(const PressureBc&)=delete;
    // Set all variables
    PressureBc(const string& segid,Pos position,const tasystem::var& p,const tasystem::var& Ts,const tasystem::var& T);
    // Assume solid temperature constant at gc.T0;
    PressureBc(const string& segid,Pos position,const tasystem::var& p,const tasystem::var& Ts); 
    // Assume above and adiabatic compresion/expansion
    PressureBc(const string& segid,Pos position,const tasystem::var& p);
    PressureBc(const PressureBc& other)=delete;

    segment::Connector* copy(const tasystem::TaSystem& s) const { return new PressureBc(*this,s);}
    virtual vd error() const;
    virtual ~PressureBc(){}
    #ifndef SWIG
    us getNEqs() const;
    virtual void updateNf();
    virtual void setEqNrs(us firsteqnr);    
    virtual void jac(tasystem::Jacobian&) const;
    // ------------------------------
    virtual void show(us i) const;
  private:
    #endif // SWIG    
  };

} // namespace tube

#endif // PRESSUREBC_H
//////////////////////////////////////////////////////////////////////

