#pragma once
#ifndef _SEGBASE_H_
#define _SEGBASE_H_
#include <vtypes.h>
#include "globalconf.h"
#include "geom.h"


namespace segment{
  class SegBase;
  using tasystem::Globalconf;
  typedef std::vector< const SegBase* > Segvec;
  
  class SegBase{
  private:
    us number=0;
    Geom* geomptr;
  public:
    Geom& geom;			// The geometry    
  protected:
    us nL=0,nR=0;
    us nleft=0,nright=0;	// Deprecated!
    string type;

  public:
    SegBase(Geom geom);
    virtual ~SegBase();
    const Globalconf* gc=NULL;	// Global configuration of the system
  protected:
    void newgeom(const Geom& newgeom);
    
  public:    
    const string& gettype() const;
    const Segvec& Right() const {return right;}
    const Segvec& Left() const {return left;}
    Segvec left,right;
    const us& getNumber() const {return number;}
    bool operator==(const SegBase& seg2) const; // Check if two segments are the same

    
    
  };
  
  
} // Namespace segment


#endif /* _SEGBASE_H_ */



