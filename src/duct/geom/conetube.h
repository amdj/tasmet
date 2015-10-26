#pragma once
#ifndef _CONETUBE_H_
#define _CONETUBE_H_

#include <exception>
#include "geom.h"
#include "vertplates.h"
#include "constants.h"


namespace duct{


 class ConeTube: public Geom{
 protected:
   d SL,SR,rL,rR;
 public:
   ConeTube(const vd& g,d r1,d r2,bool blapprox=true);
   ConeTube(const ConeTube& t):
     Geom(t),
     SL(t.SL),
     SR(t.SR),
     rL(t.rL),
     rR(t.rR)
   {}
   virtual d S(us i) const;
   virtual d phi(us i) const {return 1.0;}
   virtual d rh(us i) const;
   virtual string shape() const { return "circ";}
   virtual Geom* copy() const {return new ConeTube(*this);}
   virtual void show() const;
  };

  class CylindricalTube:public ConeTube{
    
  public:
    CylindricalTube(const vd& g,d r,bool blapprox=true);
    CylindricalTube(const CylindricalTube& t):CylindricalTube(t.x(),t.rL,t.isBlApprox()){}
    virtual Geom* copy() const {return new CylindricalTube(*this);}
    virtual void show() const;
  };


} // namespace duct


#endif /* _CONETUBE_H_ */
