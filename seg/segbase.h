#pragma once
#ifndef _SEGBASE_H_
#define _SEGBASE_H_
#include "vtypes.h"
#include "globalconf.h"
#include "geom.h"
#include "jacobian.h"


namespace segment{
  class SegBase;
  using tasystem::Globalconf;
  using tasystem::Jacobian;
  
  typedef vector<const SegBase*> SegBaseVec;

  class SegBase{
  private:
    us number=0;		// Required for TaSystem. Not used in
				// any segment code
  public:
    Geom geom;			// The geometry    
    const Globalconf* gc=NULL;	// Global configuration of the system
  protected:
    SegBaseVec left,right;
  public:
    SegBase(const Geom& geom);
    SegBase(const SegBase& o);
    SegBase& operator=(const SegBase&);
    virtual us getNDofs() const=0;
    virtual us getNEqs() const=0;    
    virtual us getNVertex() const=0;
    virtual ~SegBase(){}


    const SegBaseVec& getRight() const {return right;}
    const SegBaseVec& getLeft() const {return left;}
    void setRight(const SegBase&);	   // Couple segment to some segment on left side
    void setLeft(const SegBase&);		   // Couple segment
						   // to some segment
						   // on right side
    virtual string getType() const=0;		   // This param is
    virtual d getCurrentMass() const=0;    
    // important for connecting the segments
    
    virtual string getName() const=0; // This one is jus the name
    // ------------------------------
    virtual void init(const Globalconf&); // Implementation updates gc
    // ptr of this instance.
    virtual void setDofNrs(us firstdofnr)=0;
    virtual void setEqNrs(us firstdofnr)=0;    
    virtual vd error() const=0;
    virtual void show(us) const=0;
    virtual void jac(Jacobian&) const=0;
    virtual void domg(vd&) const=0;	// Derivative of error w.r.t. base frequency.
    virtual void dmtotdx(vd&) const=0; // Derivative of current mass in
				    // system to all dofs.
    virtual void setRes(vd res)=0;
    virtual vd getRes() const=0;
    virtual SegBase* copy() const=0;
    
    const us& getNumber() const {return number;}
    void setNumber(us number) {this->number=number;} 
    bool operator==(const SegBase& seg2) const; // Check if two segments are the same
  };
  
  
} // Namespace segment


#endif /* _SEGBASE_H_ */



