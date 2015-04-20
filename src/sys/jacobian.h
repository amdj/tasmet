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
    void operator+=(const Jacobian&);
    void operator+=(const JacRow&);
    TripletList getTriplets(us ndofs) const;
  };

  
} // namespace tasystem

#endif // JACOBIAN_H
//////////////////////////////////////////////////////////////////////
