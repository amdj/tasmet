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
#include <map>
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
  // Purge a vector of components
  template<typename Key,typename T>
  void purge(std::map<Key,T>& map){
    TRACE(10,"purge(map)");
    for (auto& it: map){
      delete it.second;
      it.second=nullptr;
    }
    map.clear();
  }
  template<typename T>
  const T& min(const T& x,const T& y)  {
    return x<=y? x : y;
  }
  template<typename T>
  const T& max(const T& x,const T& y)  {
    return x<=y? y : x;
  }
    
  
  
} // namespace utils


#endif // UTILS_H
//////////////////////////////////////////////////////////////////////
