#pragma once
#ifndef _SEGCONBASE_H_
#define _SEGCONBASE_H_
#include "vtypes.h"



namespace tasystem{
  class Jacobian;
  class Globalconf;
  class TaSystem;
}

namespace segment{
  SPOILNAMESPACE

  class SegConBase{
    static us globnr_;
    int number=-1;		// Required for TaSystem. Not used in
    // any segment/connector code
    SegConBase& operator=(const SegConBase&);
    std::string name_;
  protected:
    bool init_=false;
  public:
    const tasystem::Globalconf* gc=NULL;	// Global configuration of the system
  public:
    SegConBase();
    SegConBase(const SegConBase&);
    virtual ~SegConBase(){}
    // Get and set number
    void setNumber(us number) {this->number=number;} 
    const us& getNumber() const {return number;}

    // Get and set name
    const std::string& getName() const{return name_;} // This one is just the name
    void setName(const std::string& name){ name_=name;} // This one is just the name

    virtual bool init(const tasystem::TaSystem&); // Implementation updates gc
    bool isInit() const{return init_;}

    virtual void setEqNrs(us firstdofnr)=0;    
    virtual us getNEqs() const=0;    
    virtual vd error() const=0;
    virtual void show(us) const=0;
    virtual void jac(tasystem::Jacobian&) const=0;
    virtual void updateNf()=0;  // Update nr of frequencies

  };

} // namespace segment
#endif /* _SEGCONBASE_H_ */
