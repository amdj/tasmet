// solprogress.h
//
// Author: J.A. de Jong 
//
// Description:
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef SOLPROGRESS_H
#define SOLPROGRESS_H
#include "errorvals.h"
#include "vtypes.h"

namespace tasystem {
  
  class SolProgress
  {
    vector<ErrorVals> ev;
  public:
    SolProgress(){ ev.clear();}
    ~SolProgress(){}
    SolProgress& operator+=(const ErrorVals& e){
      ev.push_back(e);
      return *this;
    }
    vd Funer() const {
      us siz=ev.size();
      vd funer(siz);
      for(us i=0;i<siz;i++)
	funer(i)=ev.at(i).getFuner();
      return funer;
    }
  };

  
} // tasystem
#endif // SOLPROGRESS_H
//////////////////////////////////////////////////////////////////////
