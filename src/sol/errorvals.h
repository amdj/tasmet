// errorvals.h
//
// Author: J.A. de Jong 
//
// Description:
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef ERRORVALS_H
#define ERRORVALS_H
#include <atomic>
#include <cassert>
#include "vtypes.h"

namespace tasystem {
  SPOILNAMESPACE

  class ErrorVals{
    std::atomic<d> funer;
    std::atomic<d> reler;
  public:
    ErrorVals(d fe,d re): funer(fe),reler(re){}
    ErrorVals(const ErrorVals& e) {
      this->operator=(e);
    }
    d getFuner() const {return funer;}
    d getReler() const {return reler;}     
    #ifndef SWIG                // Swig does not know about
                                // initializer lists
    ErrorVals& operator=(const ErrorVals& e) {
      funer=e.getFuner();
      reler=e.getReler();
      return *this;
    }
    ErrorVals(std::initializer_list<d> il)
    {
      assert(il.size()==2);
      auto it=il.begin();
      funer=*it++;
      reler=*it;
    }
    #endif
    ~ErrorVals(){ TRACE(5,"~ErrorVals()"); }
  };

  
} // tasystem


#endif // ERRORVALS_H
//////////////////////////////////////////////////////////////////////


