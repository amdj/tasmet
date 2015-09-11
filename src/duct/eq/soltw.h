// soltw.h
//
// Author: J.A. de Jong 
//
// Description:
// Solve equation which relates the wall temperature to the
// area-averaged temperature of the solid:
//
//  Q_s->f = S_s H_s ( T_s - T_w )
//
//  where: H_s is the geometry-dependent heat transfer coefficient
//  (W/m^3K) which is computed in solidh.h.
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef SOLTW_H
#define SOLTW_H

#include "ductequation.h"
#include "constants.h"

namespace solids{
  class Solid;
}

namespace duct {
  class LaminarDuct;
  class SolidH;

  class SolTw:public Equation
  {
    const solids::Solid* solid=nullptr;
    const LaminarDuct* duct=nullptr;
    // Shape-dependent implementation of heat transfer coefficient.
  const SolidH* h=nullptr;
  public:
    SolTw(const Cell& v,const LaminarDuct& t);
    SolTw(const SolTw&)=delete;
    virtual void init();
    virtual tasystem::JacRow jac() const;
    virtual enum EqType getType() const { return EqType::SolTwEq;}
    virtual void show() const; 
    virtual vd error() const;			// Error in Energy equation at node i
    virtual void domg(vd&) const;
    virtual ~SolTw();
  };				// SolTw

} // namespace duct  


#endif // SOLTW_H
//////////////////////////////////////////////////////////////////////

