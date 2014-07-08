#pragma once
#ifndef _SEG_H_
#define _SEG_H_
#include "vertex.h"
// #include "segbase.h"
#include <vtypes.h>
#include <memory>
#include "geom.h"

#define MAXCONNECT 3		// Maximum segments that can be connected at left and right connection
namespace tube{
  class TubeVertex;
}
namespace segment{

  SPOILNAMESPACE
  typedef std::shared_ptr<Vertex> vertexptr;

  using tasystem::Globalconf;

  class Seg;
  typedef std::shared_ptr< Seg > Segptr;
  typedef std::vector< const Seg* > Segvec;
  enum SegCoupling{
    headhead,tailtail,headtail,tailhead
  };

  void coupleSegs(Seg& seg1,Seg& seg2,SegCoupling); // Couple two segments  

  class Seg{
  public:
    friend void coupleSegs(Seg&,Seg&,SegCoupling);
    friend class Vertex;
    friend class tube::TubeVertex;
    
    Seg(Geom geom); // nL,nR initiated as 0
    Seg(const Seg&);		       // Copy constructor, really copies everything
    Seg& operator=(const Seg&);
    virtual void setRight(const Seg&);	   // Couple segment to some segment on left side
    virtual void setLeft(const Seg&);		   // Couple segment to some segment on right side
    virtual void Init(const tasystem::Globalconf&);			   // Initializer method. Different for each segment type

    const Segvec& Right() const {return right;}
    const Segvec& Left() const {return left;}

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
    vd Error();			// Return error vector for this segment
    vd GetRes();		// Return result vector for this segment
    dmat Jac();			// Return Jacobian matrix
    void SetRes(vd res);
    void setnodes(us n1,us n2){ nL=n1; nR=n2;}
    virtual ~Seg(){TRACE(-5,"~Seg()");}
    

  public:
    const Globalconf* gc=NULL;	// Global configuration of the system
    Geom geom;			// The geometry
    std::vector<vertexptr> vvertex;
  protected:
    us Ndofs,Ncells;
    string type;
    us nL,nR;
    Segvec left,right;
    us nleft,nright;
  private:
    us number;

  };
  


  
} // Namespace segment


#endif /* _SEG_H_ */



