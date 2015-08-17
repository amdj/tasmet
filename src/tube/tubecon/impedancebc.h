// impedancebc.h
//
// Author: J.A. de Jong 
//
// Description:
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef IMPEDANCEBC_H
#define IMPEDANCEBC_H

#include "constants.h"
#include "tubebc.h"
#include "var.h"

#ifndef SWIG
#ifndef PyObject_HEAD
struct _object;
typedef _object PyObject;
#endif
#endif

namespace tube{
  #ifndef SWIG
  SPOILNAMESPACE
  #endif
  #ifdef SWIG
  %catches(std::exception,...) ImpedanceBc::ImpedanceBc(us segnr,Pos position,PyObject* pyfunc,d T0=constants::T0);  // %feature("notabstract") PressureBc;
  #endif // SWIG


  class ImpedanceBc:public TubeBc {
    PyObject* impedanceFunc=nullptr;
    d T0;
    ImpedanceBc(const ImpedanceBc& other,const tasystem::TaSystem&);
  public:
    // pyfunc is impedance function, which should return a complex
    // number as a function of omega: z=Z(omega)
    ImpedanceBc(us segnr,Pos pos,PyObject* pyfunc,d T0=constants::T0);
    ImpedanceBc(const ImpedanceBc&)=delete;
    virtual ~ImpedanceBc();

    segment::Connector* copy(const tasystem::TaSystem& s) const { return new ImpedanceBc(*this,s);}
    virtual vd error() const;

    #ifndef SWIG
    void updateNf();
    void setEqNrs(us firstdofnr);
    us getNEqs() const {return 3*gc->Ns();}    
    void jac(tasystem::Jacobian&) const;
    void show(us i) const;
    #endif
};



} // namespace tube


#endif /* _IMPEDANCEBC_H_ */



