#pragma once			// 
#ifndef _GEOM_H_
#define _GEOM_H_

#include "vtypes.h"
#include "localgeom.h"
#include "grid.h"
#include "pos.h"

namespace tube {
  SPOILNAMESPACE;

  class Geom{
    Grid grid_;

    Geom& operator=(const Geom& other);
    bool blapprox=false;
    bool prismatic=true;


  protected:
    Geom(const Grid& g,bool blapprox=false,bool prismatic=true);

  public:
    virtual ~Geom(){}
    LocalGeom localGeom(us i) const;	// Get a local geometry for a
                                        // certain vertex
    const Grid& grid() const{return grid_;}
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
    const vd& x() const {return grid_.getx();}
    d x(us i) const {return x()(i);}
    d L() const {return grid_.getL();}
    us gp() const {return x().size();}
    d Sleft() const {return S(0);}
    d Sright() const {return S(nCells());}
    d phileft() const {return phi(0);}
    d phiright() const {return phi(nCells());}
    d rhleft() const {return rh(0);}
    d rhright() const {return rh(nCells());}

    d vx(us i) const;		 // Vertex positions    
    vd vSf_vec() const;
    vd vS_vec() const;
    vd vphi_vec() const;
    vd vrh_vec() const;    
    vd vx_vec() const;
    d Sf(us i) const {return S(i)*phi(i);}		 // Fluid-occupied cross-sectional area of cell-wall    
    d Ss(us i) const {return (1-phi(i))*S(i);}		 // Solid-occupied cross-sectional area
    d vS(us i) const;		 // Vertex cross-sectional area
    d vSf(us i) const;		 // Fluid-occupied cross-sectional area for vertex
    d vSs(us i) const;		 // Solid-occupied cross-sectional area for vertex
    d vVf(us i) const;		 // Fluid-occupied volume of cell
    d vVs(us i) const;		 // Solid-occupied volume of cell
    
    d vphi(us i) const;			// Volume porosity of vertex
    d vrh(us i) const;			// Hydraulic radius of vertex

    d getFluidVolume() const;
    d getSolidVolume() const;


  };                            /* class Geom */
  

  
} // namespace tube

#endif /* _GEOM_H_ */
