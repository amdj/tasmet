#pragma once
#ifndef _SEG_H_
#define _SEG_H_
#include "segbase.h"
#include "vertex.h"
#include <memory>


namespace segment{

  SPOILNAMESPACE
  typedef std::unique_ptr<Vertex> vertexptr;
  class Equation;
  
  class Seg:public SegBase{
  public:
    std::vector<vertexptr> vvertex; // Owned by segment, but
				    // controlled by derived classes.
  protected:
    us nDofs=0;
  public:    
    Seg(const Geom& geom); // nL,nR initiated as 0
    Seg(const Seg&);		       // Copy constructor, really copies everything
    Seg& operator=(const Seg&);

    // Coupling of segments
    // ------------------------------ DEPRECATED!!
    
    // Initialized method (after adding to a system)
    virtual void init(const tasystem::Globalconf&);			   // Initializer method. Different for each segment type
    void cleanup();
    virtual void show(bool showVertices=false);

    const us& getNDofs() const {return nDofs;}
    const us& getNCells() const {return geom.nCells;}
    // const us& getnL() const {return nL;}
    // const us& getnR() const {return nR;}

    virtual void setLeftBc(Vertex* v); // Set left boundary condition vertex
    virtual void setRightBc(Vertex* v); // Set left boundary condition vertex    
    
    vd error() const;			// Return error vector for this segment
    vd getRes() const;		// Return result vector for this segment
    dmat jac() const;			// Return Jacobian matrix
    void setRes(vd res);
    // void setnodes(us n1,us n2){ nL=n1; nR=n2;}
    virtual ~Seg(){TRACE(-5,"~Seg()"); Seg::cleanup();}
  private:
    void showVertices();    
  };
  


  
} // Namespace segment


#endif /* _SEG_H_ */



