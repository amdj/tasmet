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

  public:
    vector<JacRow> jacrows;
    us maxRow() const;
    us maxCol() const;
    void operator+=(const Jacobian&);
    void operator+=(const JacRow&);
    TripletList getTriplets() const;
  };

  
} // namespace tasystem

#endif // JACOBIAN_H
//////////////////////////////////////////////////////////////////////
