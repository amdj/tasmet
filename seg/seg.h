#pragma once
#ifndef _SEG_H_
#define _SEG_H_
#include "segbase.h"
#include "vertex.h"
#include <memory>

#define MAXCONNECT 3		// Maximum segments that can be connected at left and right connection

namespace segment{

  SPOILNAMESPACE
  typedef std::shared_ptr<Vertex> vertexptr;

  enum SegCoupling{
    headhead,tailtail,headtail,tailhead
  };

  
  class Seg:public SegBase{
    friend void coupleSegs(Seg&,Seg&,SegCoupling);
  public:
    std::vector<vertexptr> vvertex; // Owned by segment, but controlled by derived classes.
  private:
    bool leftbc=false;
    bool rightbc=false;
  protected:
    us Ndofs=0;
    us Nvertex=0;

  public:    
    Seg(Geom geom); // nL,nR initiated as 0
    Seg(const Seg&);		       // Copy constructor, really copies everything
    Seg& operator=(const Seg&);

    // Coupling of segments
    // ------------------------------ DEPRECATED!!
    void setRight(const Seg&);	   // Couple segment to some segment on left side
    void setLeft(const Seg&);		   // Couple segment to some segment on right side
    // ------------------------------
    
    // Initialized method (after adding to a system)
    virtual void Init(const tasystem::Globalconf&);			   // Initializer method. Different for each segment type

    const us& getNdofs() const {return Ndofs;}
    const us& getNcells() const {return geom.Ncells;}
    const us& getnL() const {return nL;}
    const us& getnR() const {return nR;}

    void setLeftbc(const Vertex& v); // Set left boundary condition vertex
    void setRightbc(const Vertex& v); // Set left boundary condition vertex    
    
    vd Error();			// Return error vector for this segment
    vd GetRes();		// Return result vector for this segment
    dmat Jac();			// Return Jacobian matrix
    void SetRes(vd res);
    void setnodes(us n1,us n2){ nL=n1; nR=n2;}
    virtual ~Seg(){TRACE(-5,"~Seg()");}


  };
  


  
} // Namespace segment


#endif /* _SEG_H_ */



