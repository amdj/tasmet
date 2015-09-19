// stack.h
//
// Author: J.A. de Jong 
//
// Description:
// Laminar flow in a tube where solid heat transfer is also taken into
// account: a Stack. Can also function as a heat exchanger, when Qsin
// is given.
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef STACK_H
#define STACK_H
#include "laminarduct.h"
#include "ductwithsolid.h"

namespace duct {

  class Stack:  public DuctWithSolid, public LaminarDuct{
  public:
    Stack(const Geom& geom,const string& solid);
    Stack(const Stack&,const tasystem::TaSystem&);
    ~Stack();
    Stack& operator=(const Stack&)=delete;

    #ifndef SWIG
    segment::Seg* copy(const tasystem::TaSystem& s) const {return new Stack(*this,s);}
    void show(us) const;
    void setVarsEqs(Cell&) const;
    #endif
  };
} // namespace duct


#endif // STACK_H
//////////////////////////////////////////////////////////////////////
