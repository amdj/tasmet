#pragma once
#ifndef _CONNECTOR_H_
#define _CONNECTOR_H_
#include "vtypes.h"

namespace tasystem{
  class Jacobian;
  class Globalconf;
  class TaSystem;
}

namespace segment{
  SPOILNAMESPACE


  // A connector contains only equations, no degrees of freedom
  class Connector{
    static us globnr_;
    int number=-1;		// Required for TaSystem. Not used in
    // any segment/connector code
    Connector& operator=(const Connector&);
    std::string name_;
    bool init_=false;
  public:
    const tasystem::Globalconf* gc=NULL;	// Global configuration of the system
    Connector();
    Connector(const Connector& o);
    virtual ~Connector(){}

    // Get and set number
    void setNumber(us number) {this->number=number;} 
    const us& getNumber() const {return number;}

    // Get and set name
    const std::string& getName() const{return name_;} // This one is just the name
    void setName(const std::string& name){ name_=name;} // This one is just the name

    virtual void init(const tasystem::TaSystem&); // Implementation updates gc
    bool isInit() const{return init_;}

    virtual void setEqNrs(us firstdofnr)=0;    

    virtual us getNEqs() const=0;    
    virtual vd error() const=0;
    virtual void show(us) const=0;
    virtual void jac(tasystem::Jacobian&) const=0;
    virtual void updateNf()=0;  // Update nr of frequencies

  };
} // namespace segment

#endif /* _CONNECTOR_H_ */
