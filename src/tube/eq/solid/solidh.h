// solidh.h
//
// Author: J.A. de Jong 
//
// Description:
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef SOLIDH_H
#define SOLIDH_H
#include "vtypes.h"

namespace solids{class Solid;}
namespace tube {
  class Cell;

  class SolidH{
    SolidH* h=nullptr;
  public:
    SolidH(){}
    SolidH(const string& shape);
    virtual dmat H(const Cell& v,const solids::Solid& s) const {
      return h->H(v,s);
    }
    virtual ~SolidH(){ delete h;  }
  };
} // namespace tube

#endif // SOLIDH_H
//////////////////////////////////////////////////////////////////////
