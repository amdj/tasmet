#pragma once
#ifndef _SEG_H_
#define _SEG_H_

#include <memory>
#include "vtypes.h"


namespace tasystem{
  class Jacobian;
  class Globalconf;
  class TaSystem;
}

namespace segment{
  SPOILNAMESPACE


  class Seg{
    int number=-1;		// Required for TaSystem. Not used in
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
    tasystem::Globalconf& getGc() const;
    virtual void init(const tasystem::TaSystem&); // Implementation updates gc
    virtual Seg* copy() const=0;
    virtual void resetHarmonics()=0;
    virtual std::string getType() const=0;		   // This param is
    // important for connecting the segments
    virtual std::string getName() const=0; // This one is just the name

    // ------------------------------ config methods
    virtual void setDofNrs(us firstdofnr)=0;
    virtual void setEqNrs(us firstdofnr)=0;    
    virtual us getNDofs() const=0;
    virtual us getNEqs() const=0;    
    virtual d getCurrentMass() const=0;
    // ------------------------------
    virtual vd error() const=0;
    virtual void show(us) const=0;
    virtual void jac(tasystem::Jacobian&) const=0;
    virtual void domg(vd&) const=0;	// Derivative of error w.r.t. base frequency.
    virtual void setRes(const vd& res)=0;  // Setting result vector
    virtual void setRes(const Seg& res)=0; // Copying contents
    virtual void dmtotdx(vd&) const=0; // Derivative of current fluid mass in
    // system to all dofs.

    virtual void updateNf()=0;  // Update nr of frequencies
    virtual vd   getRes() const=0; // Get a result vector
    virtual d getRes(us dofnr) const=0; // Get result for certain dof nr
  };
  
} // Namespace segment


#endif /* _SEG_H_ */



