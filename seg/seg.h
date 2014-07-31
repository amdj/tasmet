#pragma once
#ifndef _SEG_H_
#define _SEG_H_
#include "segbase.h"
#include "vertex.h"
#include <memory>


namespace segment{

  SPOILNAMESPACE
  typedef std::unique_ptr<Vertex> vertexptr;

  
  class Seg:public SegBase{
  public:
    std::vector<vertexptr> vvertex; // Owned by segment, but
				    // controlled by derived classes.
  protected:
    us Ndofs=0;
    us Nvertex=0;
  public:    
    Seg(Geom geom); // nL,nR initiated as 0
    Seg(const Seg&);		       // Copy constructor, really copies everything
    Seg& operator=(const Seg&);

    // Coupling of segments
    // ------------------------------ DEPRECATED!!
    
    // Initialized method (after adding to a system)
    virtual void Init(const tasystem::Globalconf&);			   // Initializer method. Different for each segment type
    virtual Vertex* makeVertex(us i,const Globalconf& g);
    void show(bool showVertices=false);
    void showVertices();
    const us& getNdofs() const {return Ndofs;}
    const us& getNcells() const {return geom.Ncells;}
    // const us& getnL() const {return nL;}
    // const us& getnR() const {return nR;}

    void setLeftbc(Vertex* v); // Set left boundary condition vertex
    void setRightbc(Vertex* v); // Set left boundary condition vertex    
    
    vd Error();			// Return error vector for this segment
    vd GetRes();		// Return result vector for this segment
    vd GetResAt(us varnr,us freqnr); // Extract a result vector for given variable number (rho,U,T,p,Ts) and frequency number.
    dmat Jac();			// Return Jacobian matrix
    void SetRes(vd res);
    // void setnodes(us n1,us n2){ nL=n1; nR=n2;}
    virtual ~Seg(){TRACE(-5,"~Seg()");}

  };
  


  
} // Namespace segment


#endif /* _SEG_H_ */



