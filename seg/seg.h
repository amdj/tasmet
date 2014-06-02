#pragma once
#ifndef _SEG_H_
#define _SEG_H_
#include "vertex.h"
#include <vtypes.h>
#include <memory>

namespace segment{

  SPOILNAMESPACE
  typedef std::shared_ptr<Vertex> vertexptr;
  
  class Seg;

  enum SegCoupling{
    headhead,tailtail,headtail,tailhead
  };

  void coupleSegs(Seg& seg1,Seg& seg2,SegCoupling); // Couple two segments  
  
  class Seg{
  public:
    friend void coupleSegs(Seg&,Seg&,SegCoupling);
    
    Seg(const tasystem::Globalconf& gc); // nL,nR initiated as 0
    virtual void setRight(const Seg&);	   // Couple segment to some segment on left side
    virtual void setLeft(const Seg&);		   // Couple segment to some segment on right side
    virtual void Init()=0;			   // Initializer method. Different for each segment type
    const Seg* Right() const {return right;}
    const Seg* Left()const {return left;}
    const us& getNumber() const {return number;}
    const us& getNdofs() const {return Ndofs;}
    const us& getNcells() const {return Ncells;}
    const us& getnL() const {return nL;}
    const us& getnR() const {return nR;}


    void setLeftbc(vertexptr v); // Set left boundary condition vertex
    void setRightbc(vertexptr v); // Set left boundary condition vertex    
    void setLeftbc(Vertex* v); // Set left boundary condition vertex, takes over ownership of object
    void setRightbc(Vertex* v); // Set right boundary condition vertex, takes over ownership of object
    
    
    const string& gettype() const {return type;}
    bool operator==(const Seg& seg2) const; // Check if two segments are the same

    const tasystem::Globalconf& gc;	// Global configuration of the system
    vd Error();			// Return error vector for this segment
    vd GetRes();		// Return result vector for this segment
    dmat Jac();			// Return Jacobian matrix
    void SetRes(vd res);
    void setnodes(us n1,us n2){ nL=n1; nR=n2;}

    const us& Ns;
    virtual ~Seg() {}

  protected:
    std::vector<vertexptr> vvertex;
    us Ndofs,Ncells;
    string type;
    us nL,nR;
    Seg const *left,*right;    
  private:
    us number;
  };


  
} // Namespace segment


#endif /* _SEG_H_ */



