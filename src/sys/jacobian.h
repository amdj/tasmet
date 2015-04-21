// jacobian.h
//
// Author: J.A. de Jong 
//
// Description:
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef JACOBIAN_H
#define JACOBIAN_H

#include "jacrow.h"

namespace tasystem{
  SPOILNAMESPACE
  class TripletList;

  class Jacobian{
    us ndofs_;
  public:
    Jacobian(us ndofs):ndofs_(ndofs){}
    vector<JacRow> jacrows;
    void operator+=(const Jacobian&);
    void operator+=(const JacRow&);
    operator TripletList() const;
  };

  
} // namespace tasystem

#endif // JACOBIAN_H
//////////////////////////////////////////////////////////////////////
