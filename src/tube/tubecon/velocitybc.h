// velocitybc.h
//
// Author: J.A. de Jong 
//
// Description:
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef VELOCITYBC_H
#define VELOCITYBC_H

#include "constants.h"
#include "tubebc.h"
#include "var.h"

namespace tube{
  #ifndef SWIG
  SPOILNAMESPACE
  #endif
  #ifdef SWIG
  %catches(std::exception,...) VelocityBc::VelocityBc(const string& segid,Pos position,const tasystem::var& u,\
                                                      d T0=constants::T0,bool arbitrateMass=false);  // %feature("notabstract") PressureBc;
  #endif // SWIG


  class VelocityBc:public TubeBc {
    bool arbitrateMass=false;
    d T0;                       // Reference temperature for adiabatic
                                // compression/expansion 
    tasystem::var u_p;          // Prescribed velocity
    VelocityBc(const VelocityBc& other,const tasystem::TaSystem&);
  public:
    // u is prescribed velocity field
    VelocityBc(const string& segid,Pos pos,const tasystem::var& u,d T0=constants::T0,bool arbitrateMass=false);
    VelocityBc(const VelocityBc&)=delete;
    virtual ~VelocityBc();

    segment::Connector* copy(const tasystem::TaSystem& s) const { return new VelocityBc(*this,s);}
    virtual vd error() const;
    int arbitrateMassEq() const;
    #ifndef SWIG
    void updateNf();
    void setEqNrs(us firstdofnr);
    us getNEqs() const {return 3*gc->Ns();}    
    void jac(tasystem::Jacobian&) const;
    void show(us i) const;
    #endif
};



} // namespace tube

#endif // VELOCITYBC_H
//////////////////////////////////////////////////////////////////////
