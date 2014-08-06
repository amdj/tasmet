#pragma once
#ifndef _SEGBASE_H_
#define _SEGBASE_H_
#include <vtypes.h>
#include "globalconf.h"
#include "geom.h"

namespace segment{
  class SegBase;
  using tasystem::Globalconf;
  typedef vector<const SegBase*> SegBaseVec;

  class SegBase{
  private:
    us number=0;
  public:
    Geom geom;			// The geometry    
    const Globalconf* gc=NULL;	// Global configuration of the system
  protected:
    SegBaseVec left,right;
    string type;
    
  public:
    SegBase(const Geom& geom);
    virtual ~SegBase(){ cleanup();};
    void cleanup(){}		   // Stub method in case class contains any dynamic allocated data
    void setRight(const SegBase&);	   // Couple segment to some segment on left side
    void setLeft(const SegBase&);		   // Couple segment to some segment on right side
    // ------------------------------
    virtual void init(const Globalconf&);
    const us& getNumber() const {return number;}
    void setNumber(us number) {this->number=number;} 
    const string& gettype() const;
    const SegBaseVec& getRight() const {return right;}
    const SegBaseVec& getLeft() const {return left;}
    bool operator==(const SegBase& seg2) const; // Check if two segments are the same
  };
  
  
} // Namespace segment


#endif /* _SEGBASE_H_ */



