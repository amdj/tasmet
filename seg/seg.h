#pragma once
#ifndef _SEG_H_
#define _SEG_H_
#include "segbase.h"
#include <memory>


namespace segment{
  SPOILNAMESPACE

  using std::string;

  namespace geom{
    class Geom;
  }
  namespace tasystem{
    class Jacobian;
    class Globalconf;
  }

  class Seg{
  private:
    us number=0;		// Required for TaSystem. Not used in
    // any segment code
    Seg& operator=(const Seg&);
  public:
    const Globalconf* gc=NULL;	// Global configuration of the system
  public:
    Seg(){}
    Seg(const Seg& o);
    const us& getNumber() const {return number;}
    void setNumber(us number) {this->number=number;} 
    bool operator==(const Seg& other) const {return (this==&other);}
    virtual Seg* copy() const=0;
    virtual void resetHarmonics()=0;
    virtual void init(const Globalconf&); // Implementation updates gc
    virtual us getNDofs() const=0;
    virtual us getNEqs() const=0;    
    virtual ~Seg(){}            // We do not own the gc instance
    const Geom& geom() const {return *geom_;}
    virtual string getType() const=0;		   // This param is
    virtual d getCurrentMass() const=0;    
    // important for connecting the segments
    virtual string getName() const=0; // This one is jus the name
    // ------------------------------
    // ptr of this instance.
    virtual void setDofNrs(us firstdofnr)=0;
    virtual void setEqNrs(us firstdofnr)=0;    
    virtual vd error() const=0;
    virtual void show(us) const=0;
    virtual void jac(Jacobian&) const=0;
    virtual void domg(vd&) const=0;	// Derivative of error w.r.t. base frequency.
    virtual void dmtotdx(vd&) const=0; // Derivative of current fluid mass in
    // system to all dofs.
    virtual void setRes(vd res)=0;
    virtual void updateNf()=0;
    virtual void setRes(const Seg&)=0;    
    virtual vd   getRes() const=0;
  };
  
} // Namespace segment


#endif /* _SEG_H_ */



