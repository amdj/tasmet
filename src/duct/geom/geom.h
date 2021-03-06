#pragma once			// 
#ifndef _GEOM_H_
#define _GEOM_H_

#include "constants.h"
#include "vtypes.h"
#include "localgeom.h"


namespace duct {
  #ifndef SWIG
  SPOILNAMESPACE;
  #endif

  class Geom{
    Geom& operator=(const Geom& other);
    vd x_;			// The grid
    bool blapprox=false;
    bool prismatic=true;
  protected:
    Geom(const vd& x,bool blapprox=false,bool prismatic=true);
    Geom(const Geom&);
  public:
    virtual ~Geom(){}
    #ifndef SWIG
    LocalGeom localGeom(us i) const;	// Get a local geometry for a
                                        // certain cell
    #endif
    virtual void show() const=0;
    virtual d S(us i) const=0;		 // Cross sectional area as a function of x
    virtual d phi(us i) const=0;		 // Volume porosity at position of cell walls
    virtual d rh(us i) const=0;		 // Hydraulic radius
    virtual Geom* copy() const=0;
    // Shape keyword: currently available: 'circ','vert'
    virtual string shape() const=0;

    void setPrismatic(bool isprismatic){prismatic=isprismatic;}
    bool isPrismatic() const {return prismatic;}
    void setBlApprox(bool b){ blapprox=b;}
    bool isBlApprox() const {return blapprox;}

    us nCells() const {return x().size()-1;}
    const vd& x() const {return x_;}
    d x(us i) const {return x()(i);}
    d L() const {return *(x_.end()-1);}
    us gp() const {return x().size();}
    d Sleft() const {return S(0);}
    d Sright() const {return S(nCells());}
    d phileft() const {return phi(0);}
    d phiright() const {return phi(nCells());}
    d rhleft() const {return rh(0);}
    d rhright() const {return rh(nCells());}

    d vx(us i) const;		 // Cell positions    
    vd vSf_vec() const;
    vd vS_vec() const;
    vd vphi_vec() const;
    vd vrh_vec() const;    
    vd vx_vec() const;
    d Sf(us i) const {return S(i)*phi(i);}		 // Fluid-occupied cross-sectional area of cell-wall    
    d Ss(us i) const {return (1-phi(i))*S(i);}		 // Solid-occupied cross-sectional area
    d vS(us i) const;		 // Cell cross-sectional area
    d vSf(us i) const;		 // Fluid-occupied cross-sectional area for cell
    d vSs(us i) const;		 // Solid-occupied cross-sectional area for cell
    d vVf(us i) const;		 // Fluid-occupied volume of cell
    d vVs(us i) const;		 // Solid-occupied volume of cell
    
    d vphi(us i) const;			// Volume porosity of cell
    d vrh(us i) const;			// Hydraulic radius of cell

    d getFluidVolume() const;
    d getSolidVolume() const;


  };                            /* class Geom */
  

  
} // namespace duct

#endif /* _GEOM_H_ */
