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
namespace duct {
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
} // namespace duct

#endif // SOLIDH_H
//////////////////////////////////////////////////////////////////////
