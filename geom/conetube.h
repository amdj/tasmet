#pragma once
#ifndef _CONETUBE_H_
#define _CONETUBE_H_

#include "geom.h"

namespace geom{


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

} // namespace geom


#endif /* _CONETUBE_H_ */
