#pragma once
#ifndef _SEG_H_
#define _SEG_H_

#include <memory>
#include "vtypes.h"


namespace tasystem{
  class Jacobian;
  class Globalconf;
}

namespace segment{
  SPOILNAMESPACE

  using std::string;

  class Seg{
    us number=0;		// Required for TaSystem. Not used in
    // any segment code
    Seg& operator=(const Seg&);
    bool init_=false;
  public:
    const tasystem::Globalconf* gc=NULL;	// Global configuration of the system

    Seg();
    Seg(const Seg& o);
    virtual ~Seg(){}            // We do not own the gc instance
    const us& getNumber() const {return number;}
    void setNumber(us number) {this->number=number;} 
    bool operator==(const Seg& other) const {return (this==&other);}
    bool isInit() const{return init_;}
    // Pure virtual functions
    virtual void init(const tasystem::Globalconf&); // Implementation updates gc
    virtual Seg* copy() const=0;
    virtual void resetHarmonics()=0;
    virtual string getType() const=0;		   // This param is
    virtual d getCurrentMass() const=0;    
    // important for connecting the segments
    virtual string getName() const=0; // This one is just the name
    // ------------------------------ config methods
    virtual void setDofNrs(us firstdofnr)=0;
    virtual void setEqNrs(us firstdofnr)=0;    
    virtual us getNDofs() const=0;
    virtual us getNEqs() const=0;    
    // ------------------------------
    virtual vd error() const=0;
    virtual void show(us) const=0;
    virtual void jac(tasystem::Jacobian&) const=0;
    virtual void domg(vd&) const=0;	// Derivative of error w.r.t. base frequency.
    virtual void setRes(vd res)=0;
    virtual void dmtotdx(vd&) const=0; // Derivative of current fluid mass in
    // system to all dofs.

    virtual void updateNf()=0;
    virtual void setRes(const Seg&)=0;    
    virtual vd   getRes() const=0;
    virtual d getRes(us dofnr) const=0;
  };
  
} // Namespace segment


#endif /* _SEG_H_ */



