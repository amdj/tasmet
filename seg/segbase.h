#pragma once
#ifndef _SEGBASE_H_
#define _SEGBASE_H_
#include <vtypes.h>
#include "globalconf.h"
#include "geom.h"
#define Neq (5)

namespace segment{
  class SegBase;
  using tasystem::Globalconf;
  typedef vector<const SegBase*> SegBaseVec;

  class SegBase{
  private:
    us number=0;		// Required for TAsystem. Not used in
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
    virtual ~SegBase(){ cleanup();}
    void cleanup(){}		   // Stub method in case class contains any dynamic allocated data
    void setRight(const SegBase&);	   // Couple segment to some segment on left side
    void setLeft(const SegBase&);		   // Couple segment
						   // to some segment
						   // on right side
    virtual string getType() const=0;
    // ------------------------------
    virtual void init(const Globalconf&); // Implementation updates gc
    // ptr of this instance.
    virtual vd error() const=0;
    virtual void show(bool) const=0;
    virtual dmat jac() const=0;
    virtual void setRes(vd res)=0;
    virtual vd getRes() const=0;
    virtual SegBase* copy() const=0;
    
    const us& getNumber() const {return number;}
    void setNumber(us number) {this->number=number;} 
    const SegBaseVec& getRight() const {return right;}
    const SegBaseVec& getLeft() const {return left;}
    bool operator==(const SegBase& seg2) const; // Check if two segments are the same
  };
  
  
} // Namespace segment


#endif /* _SEGBASE_H_ */



