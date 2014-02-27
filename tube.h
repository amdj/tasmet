/*
 * tube.h
 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */
#pragma once

#ifndef TUBE_H_
#define TUBE_H_
#include "globalconf.h"
#include "var/var.h"
#include "common/vtypes.h"
#include "common/material.h"


namespace segment{

  class Seg{
  public:

    Seg(globalconf::Globalconf& gc); // nL,nR initiated as 0
    ~Seg(){}
    us nL,nR;
    globalconf::Globalconf& gc;	// Global configuration of the system
    void setnodes(us n1,us n2){ nL=n1; nR=n2;}
  };

}

namespace tube {
  class Geom{
  public:
    Geom(us gp,d L,d S,d phi,d rh,string cshape);
    ~Geom();

    const us gp;		 /* Number of gridpoints */
    d L;			 /* Length of the Tube */
    vd x;			/* Grid points vector */
    vd S;			/* Cross sectional area as a function of x */
    vd Sf;			/* Fluid-occupied cross-sectional area */
    vd Ss;			/* Solid-occupied cross-sectional area */
    vd phi;			/* Volume porosity */
    vd rh;			/* Hydraulic radius */
    string shape;		/* Cross-sectional area shape keyword */

  };				/* class Geom */
  // class SubJacobian{ //Create a square matrix with values for the SubJacobian of equation x to variable y
  //   // For example, the derivative of the continuity equation at node i to the variable rho at node i+1
  // public:
  //   SubJacobian(TubeGp& tubegp);
  //   ~SubJacobian();
  //   void set(dmat);
  //   TubeGp& tubegp;		// Reference to current gridpoint
  //   us node;			// Node number
  //   dmat SubJac;		// The submatrix jacobian
  // };
  class TubeGp;
  class TubeLc{	// Tube continuity equation 
  public:
    TubeLc(TubeGp& gp);
    ~TubeLc();
    us i; 			// Current node
    d Sf,dxp,dxm;		// Fluid cs-area, dx+,dx-
    dmat operator()();

    TubeGp& tubegp;		// Reference to parent (current gridpoint)
    dmat rhoip1;	// Derivative of continuity equation to density at node i + 1
    dmat rhoi;	// Derivative of continuity equation to density at node i
    dmat rhoim1;	// Derivative of continuity equation to density at node i - 1
    dmat Uip1;	// Derivative of continuity equation to Volume flow at node i + 1
    dmat Uim1;	// Derivative of continuity equation to Volume flow at node i - 1


};				// TubeLc class

  class Tube;
  class TubeGp{ //Gridpoint at a position in a Tube
  public:
    TubeGp(const Tube& tube,us i);
    ~TubeGp();
    const Tube& tube;			// Reference to parent tube
    us i;			// The node number of this gridpoint
    TubeLc lc;			// Continuity equation
    // TubeLm lm; 			// Momentum equation

    variable::var rho;		// Density
    variable::var U;		// Volume flow
    variable::var T;		// Temperature
    variable::var p;		// Pressure
  };				// TubeGp class

  class Tube:public segment::Seg {
  public:
    Tube(globalconf::Globalconf& g,Geom geom);
    ~Tube();
    void Init(d T0,d p0);
    variable::varoperations vop;
    Geom geom;			// The geometry
    gases::Gas& gas;		// The gas in the system. Reference variable to globalconf.gas
    vector<TubeGp> gps;		// Vector of gridpoints

  protected:

    friend class TubeGp;
    friend class TubeLc;

  };				// Tube class


} /* namespace tube */



#endif /* TUBE_H_ */
