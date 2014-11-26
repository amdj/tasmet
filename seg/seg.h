#pragma once
#ifndef _SEG_H_
#define _SEG_H_
#include "segbase.h"
#include <memory>


namespace segment{
  SPOILNAMESPACE
  class Vertex;  
  class Seg:public SegBase{
  public:
    // This vector is common in each derived segment. It is controlled
    // by the derived segment as well. The segment itself does do
    // nothing with it.
    std::vector<Vertex*> vvertex;

  public:    
    Seg(const Geom& geom);
    Seg(const Seg&);		       // Copy constructor, really copies everything
    Seg& operator=(const Seg&);

    // Coupling of segments
    // Initialized method (after adding to a system)
    virtual void init(const tasystem::Globalconf&);			   // Initializer method. Different for each segment type
    void cleanup();
    void show(us showVertices=0) const;
    virtual us getNVertex() const {return vvertex.size();}    
    void domg(vd&) const;
    vd error() const;			// Return error vector for this segment
    vd getRes() const;		// Return result vector for this segment
    void jac(Jacobian&) const;			// Return Jacobian matrix
    void setRes(vd res);
    // void setnodes(us n1,us n2){ nL=n1; nR=n2;}
    virtual ~Seg(){TRACE(-5,"~Seg()");}
  private:
    void showVertices() const ;    
  };
  


  
} // Namespace segment


#endif /* _SEG_H_ */



