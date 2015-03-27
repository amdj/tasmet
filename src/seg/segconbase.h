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

  class SegConBase{
    static us globnr_;
    int number=-1;		// Required for TaSystem. Not used in
    // any segment/connector code
    std::string name_;
    bool init_=false;
  public:
    const tasystem::Globalconf* gc=NULL;	// Global configuration of the system
  public:
    SegConBase();
    SegConBase(const SegConBase&);
    SegConBase& operator=(const SegConBase&)=delete;
    virtual ~SegConBase(){}

    void setInit(bool init){init_=init;}
    bool isInit() const{return init_;}
    void checkInit() const throw(std::exception) {if(!isInit()) throw MyError("Not initialized"); }
    virtual vd error() const=0;

    // Get and set name
    const std::string& getName() const{return name_;} // This one is just the name
    void setName(const std::string& name){ name_=name;} // This one is just the name

    #ifndef SWIG
    // Get and set number
    void setNumber(us number) {this->number=number;} 
    const int& getNumber() const {return number;}


    virtual bool init(const tasystem::TaSystem&); // Implementation updates gc

    virtual void setEqNrs(us firstdofnr)=0;    
    virtual us getNEqs() const=0;    
    virtual void show(us) const=0;
    virtual void jac(tasystem::Jacobian&) const=0;
    virtual void updateNf()=0;  // Update nr of frequencies
    #endif
  };

} // namespace segment
#endif /* _SEGCONBASE_H_ */
