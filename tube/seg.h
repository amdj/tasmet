#pragma once
#ifndef _SEG_H_
#define _SEG_H_
#include "vertex.h"
#include <vtypes.h>
#include <memory>

namespace segment{
typedef Eigen::SparseMatrix<double> SpMat;
  SPOILNAMESPACE
  typedef std::shared_ptr<Vertex> vertexptr;
  
  class Seg;

  enum SegCoupling{
    headhead,tailtail,headtail,tailhead
  };

  void coupleSegs(Seg* seg1,Seg* seg2,SegCoupling); // Couple two segments
  
  class Seg{
  public:
    friend void coupleSegs(Seg*,Seg*,SegCoupling);
    
    Seg(tasystem::Globalconf& gc); // nL,nR initiated as 0
    void setRight(Seg*);	   // Couple segment to some segment on left side
    void setLeft(Seg*);		   // Couple segment to some segment on right side
    bool operator==(const Seg& seg2); // Check if two segments are the same
    us nL,nR;
    tasystem::Globalconf& gc;	// Global configuration of the system
    vd Error();
    vd GetRes();
    dmat Jac();		// Sparse matrix
    void SetRes(vd res);
    void setnodes(us n1,us n2){ nL=n1; nR=n2;}
    std::vector<vertexptr> vvertex;
    // Vertex** vvertex; // Vector of vertices
    us Ndofs,Ncells;
    const us& Ns;
  protected:
    Seg *left,*right;    
  private:
    us number;
  };

  
  
} // Namespace segment


#endif /* _SEG_H_ */
