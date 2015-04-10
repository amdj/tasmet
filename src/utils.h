// utils.h
//
// Author: J.A. de Jong 
//
// Description:
// Some generic utils.
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef UTILS_H
#define UTILS_H
#include <vector>
#include <map>
#include "tracer.h"
#include <typeinfo>
#include <exception>

namespace utils {
  SPOILNAMESPACE
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
  template<typename SegType,typename Sys>
  SegType* copySeg(const SegType& t,const Sys& sys) {
    SegType* newt=new SegType(t);
    if(!newt){
      WARN("Copying " << typeid(t).name() << "failed!");
    }
    return newt;
  } // copySeg
  
  
} // namespace utils


#endif // UTILS_H
//////////////////////////////////////////////////////////////////////

