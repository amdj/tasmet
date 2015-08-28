// mechbc.h
//
// Author: J.A. de Jong 
//
// Description:
// Mechanical domain boundary condition for a piston. With this 'connector',
// a mechanical boundary condition can be set for the piston force, or
// the piston displacement. Three types of relations can be set:
// 1: A force boundary condition
// 2: A prescribed piston displacement
// 3: A prescribed mechanical impedance
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef MECHBC_H
#define MECHBC_H
#include "connector.h"
#include "var.h"
#include "constants.h"

namespace mech {
  class Piston;
  #ifdef SWIG
  %catches(std::exception,...) MechBc::MechBc(const string& segid,Varnr var,const tasystem::var& bc);
  #endif

  class MechBc:public segment::Connector{
    // The segid should be a Piston!
    string segid;
    const Piston *p=nullptr;
    // Type of b.c. (F,x, or Z)
    Varnr type;
    // In which the boundary condition is stored:
    tasystem::var bc;
    us firsteqnr;
  public:
    MechBc(const string& segid,Varnr var,const tasystem::var& bc);
    MechBc(const MechBc&,const tasystem::TaSystem&);
    ~MechBc();
    vd error() const;
    void jac(tasystem::Jacobian&) const;
    segment::Connector* copy(const tasystem::TaSystem& s) const{return new MechBc(*this,s);}

    // Number the internal equations
    void setEqNrs(us firsteqnr){this->firsteqnr=firsteqnr;}
    // Return the total number of equations in this segment/connector.
    us getNEqs() const;    
    void show(us) const;
    void updateNf();  // Update nr of frequencies.

  };
  
} // namespace mech

#endif // MECHBC_H
//////////////////////////////////////////////////////////////////////
