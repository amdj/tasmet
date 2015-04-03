// utils.h
//
// Author: J.A. de Jong 
//
// Description:
// Some generic utils which safes typing code
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef UTILS_H
#define UTILS_H
#include <vector>
#include "tracer.h"

namespace utils {
  
  // Purge a vector of components
  template<typename T>
  void purge(std::vector<T>& vec){
    TRACE(10,"purge(vector)");
    for (T& it: vec){
      delete it;
      it=nullptr;
    }
    vec.clear();
  }
  
  
  
} // namespace utils


#endif // UTILS_H
//////////////////////////////////////////////////////////////////////
