// ductwithsolid.h
//
// Author: J.A. de Jong 
//
// Description:
// When a piece of duct also includes a solid, this class has to be
// used. The class is derived virtual public, which means that in
// derives on Duct. Hower, the Duct constructor is not called by this
// class.
// 
// WARNING: This class expects that setVarsEqs is called AFTER
// setVarsEqs of Duct is called. Moreover, Duct::setVarsEqs is not
// called by DuctWithSolid::setVarsEqs.
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef DUCTWITHSOLID_H
#define DUCTWITHSOLID_H
#include "duct.h"

namespace solids{
  class Solid;
}
  #ifdef SWIG
  %catches(std::exception,...) DuctWithSolid::DuctWithSolid(const string& solid);
  %catches(std::exception,...) DuctWithSolid::setQsin(d Qsin);
  #endif

namespace duct {
  class DuctWithSolid: virtual public Duct{
    const solids::Solid* solid=nullptr;
    d Qsin=0;			// Heat input into the solid
    d ksfrac;			// Input foezelfactor for conduction
  public:
    DuctWithSolid(const string& solid,d ksfrac=1.0);
    DuctWithSolid(const DuctWithSolid&);
    ~DuctWithSolid();
    // This one overrides Duct's hasSolid(), which would return false
    bool hasSolid() const {return true;}
    const solids::Solid& getSolid() const { return *solid;}
    // Apply solid heat input. For this, a
    // solid need to be defined with
    void setQsin(d Qsin);
    void setVarsEqs(Cell& c) const;
  };
} // namespace duct


#endif // DUCTWITHSOLID_H
//////////////////////////////////////////////////////////////////////
