/*
 * tube.h
 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#pragma once
#ifndef TUBE_H_
#define TUBE_H_
#include "../globalconf.h"
#include "geom.h"
#include "drag.h"
#include "vertex.h"
#include "../var/var.h"
#include <vtypes.h>
#include <material.h>
#include <math_common.h>

#include <Eigen/SparseCore>



namespace segment{
typedef Eigen::SparseMatrix<double> SpMat;

  class Seg{
  public:
    Seg(tasystem::Globalconf& gc); // nL,nR initiated as 0
    virtual ~Seg();
    us nL,nR;
    tasystem::Globalconf& gc;	// Global configuration of the system
    variable::varoperations vop;
    vd Error();
    vd GetRes();
    sdmat Jac();		// Sparse matrix
    void SetRes(vd res);
    
    void setnodes(us n1,us n2){ nL=n1; nR=n2;}
    Vertex** vvertex; // Vector of vertices
    us Ndofs,Ncells;
  };

} // Namespace segment

namespace tube{

  using arma::sp_mat;  
  using segment::Seg;
  using segment::Vertex;
  class Tube:public Seg {
  public:
    Tube(tasystem::Globalconf& g,Geom geom);
    Tube(const Tube& othertube); // Copy constructor to copy vertex
				 // vector ofpointers
    ~Tube();

    void setLeftbc(TubeVertex* v); // Set left boundary condition vertex
    void setRightbc(TubeVertex* v); // Set left boundary condition vertex    
    Geom geom;			// The geometry
    gases::Gas& gas;		// The gas in the system. Reference variable to globalconf.gas
    void DoIter(d dampfac=1.0);		// And damp with a factor
    vd GetResAt(us varnr,us freqnr);

  protected:
    void Init();
    LaminarDragResistance drag;
    friend class TubeVertex;
    friend class Continuity;
    friend class Momentum;
  };				// Tube class

} /* namespace tube */

#endif /* TUBE_H_ */






