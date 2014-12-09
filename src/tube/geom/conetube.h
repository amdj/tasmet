#pragma once
#ifndef _CONETUBE_H_
#define _CONETUBE_H_

#include "geom.h"
#include "vertplates.h"
#include "pos.h"
namespace tube{


 class ConeTube: public Geom{
 protected:
   d SL,SR,rL,rR;
 public:
   ConeTube(const Grid& g,d r1,d r2,bool blapprox=true);
   ConeTube(const ConeTube& t): ConeTube(t.grid(),t.rL,t.rR,t.isBlApprox()){}
   virtual d S(us i) const;
   virtual d phi(us i) const {return 1.0;}
   virtual d rh(us i) const;
   virtual string shape() const { return "circ";}
   virtual Geom* copy() const {return new ConeTube(*this);}
   virtual void show() const;
  };

  class CylindricalTube:public ConeTube{
    
  public:
    CylindricalTube(const Grid& g,d r,bool blapprox=true);
    CylindricalTube(const CylindricalTube& t):CylindricalTube(t.grid(),t.rL,t.isBlApprox()){}
    virtual Geom* copy() const {return new CylindricalTube(*this);}
    virtual void show() const;
  };

  

  class TransitionCylindricalTube:public Geom{
    Transition transition;
    d r_,perc,S_;
    d S_other,phi_other,rh_other;
  public:
    TransitionCylindricalTube(const Grid& g,d r,\
                              segment::pos TransitionSide,const Geom& other,\
                              segment::pos sideofremote,d perc=10,bool blapprox=true);
    // Copy constructor only needs to copy all params, so default will
    // suffice.
    virtual ~TransitionCylindricalTube(){}
    virtual d S(us i) const;
    virtual d phi(us i) const;
    virtual d rh(us i ) const;
    virtual string shape() const { return "circ";}
    virtual Geom* copy() const {return new TransitionCylindricalTube(*this);}
    virtual void show() const;

  };


} // namespace tube


#endif /* _CONETUBE_H_ */
