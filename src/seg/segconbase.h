// segconbase.h
//
// Author: J.A. de Jong 
//
// Description: The SegConBase class is the base class for both the
// segments and the connectors. The main difference between a segment
// and a connector is, that a connector only provide equations, but no
// degrees of freedom. A segment contains both. All common stuff, such as a name, initialization 
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef _SEGCONBASE_H_
#define _SEGCONBASE_H_
#include "vtypes.h"
#include "exception.h"
#include "phaseconstraint.h"

namespace tasystem{
  class Jacobian;
  class Globalconf;
  class TaSystem;
}

namespace segment{
  #ifndef SWIG
  SPOILNAMESPACE
  #endif

  #ifdef SWIG
  %catches(std::exception,...) SegConBase::setName();
  %catches(std::exception,...) SegConBase::init();
  #endif // SWIG

  class SegConBase{
    static us globnr_;
    int number=-1;		// Required for TaSystem. Not used in
    // any segment/connector code
    std::string name_;
    std::string id_;
  protected:
    const tasystem::Globalconf* gc=nullptr;	// Global configuration of the system
    SegConBase();
    SegConBase(const SegConBase&,const tasystem::TaSystem&);
    SegConBase(const SegConBase&)=delete;
    SegConBase& operator=(const SegConBase&)=delete;
  public:
    virtual ~SegConBase(){}

    const char* __repr__() const {return getName().c_str();}

    // Get and set name
    const std::string& getName() const{return name_;} // This one is just the name
    const std::string& getID() const{return id_;} // This one is just the name
    void setName(const std::string& name){ name_=name;} // This one is just the name
    void setID(const std::string& id){ id_=id;} // Set ID
    
    // Tell a TaSystem whether this Segment of Connector arbitrates
    // Mass or not. The special return value of -1 tells it does
    // not. If it does, the derived class should return which equation
    // should be overwritten with the mass arbitration equation.
    virtual int arbitrateMassEq() const {return -1;}
    virtual vd error() const=0;
    #ifndef SWIG
    // Number the internal equations
    virtual void setEqNrs(us firstdofnr)=0;    
    // Return the total number of equations in this segment/connector.
    virtual us getNEqs() const=0;    
    virtual void show(us) const=0;
    // Fill Jacobian with values from the equations in this
    // segment/connector.
    virtual void jac(tasystem::Jacobian&) const=0;
    virtual void updateNf()=0;  // Update nr of frequencies.

    // Return reference to Glocalconf
    const tasystem::Globalconf& Gc() const {return *gc;}

    // After the constructor is called, the initialization function is
    // called by the TaSystem object. This is done because sometimes
    // the initialization requires calling virtual functions. When
    // something goes wrong, the function should throw an
    // exception. When this happens, the TaSystem class will delete
    // the allocated segment or connector and rethrow the error.
    virtual void init() {}

    // Get and set number. These functions are accessed from TaSystem
    void setNumber(us number);
    const int& getNumber() const {return number;}



    #endif
  };

} // namespace segment
#endif /* _SEGCONBASE_H_ */
//////////////////////////////////////////////////////////////////////
