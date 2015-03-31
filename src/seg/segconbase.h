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
  #endif // SWIG

  class SegConBase{
    static us globnr_;
    int number=-1;		// Required for TaSystem. Not used in
    // any segment/connector code
    std::string name_;
    bool init_=false;
  protected:
    const tasystem::Globalconf* gc=nullptr;	// Global configuration of the system
    SegConBase();
    SegConBase(const SegConBase&);
    SegConBase& operator=(const SegConBase&)=delete;
  public:
    virtual ~SegConBase(){}

    // Get and set name
    const std::string& getName() const{return name_;} // This one is just the name
    void setName(const std::string& name){ name_=name;} // This one is just the name

    // This function determines whether a class is abstract or not for SWIG
    // Return error from internal equations
    virtual string getType() const=0;
    #ifndef SWIG
    virtual vd error() const=0;
    virtual void setEqNrs(us firstdofnr)=0;    
    virtual us getNEqs() const=0;    
    virtual void show(us) const=0;
    virtual void jac(tasystem::Jacobian&) const=0;
    virtual void updateNf()=0;  // Update nr of frequencies

    // Return reference to Glocalconf
    const tasystem::Globalconf& Gc() const {return *gc;}

    // Initialization functions

    // Initialize the Segment or connecter. Should only be done once
    // it is part of a TaSystem. Not doing this will result in
    // SegFaults, as pointers to this TaSystem will be saved.
    virtual void init(const tasystem::TaSystem&); // Implementation updates gc
    void setInit(bool init){init_=init;}
    bool isInit() const{return init_;}
    void checkInit() const {if(!isInit()) throw MyError("Not initialized"); }

    // Get and set number. These functions are accessed from TaSystem
    void setNumber(us number);
    const int& getNumber() const {return number;}



    #endif
  };

} // namespace segment
#endif /* _SEGCONBASE_H_ */
//////////////////////////////////////////////////////////////////////
