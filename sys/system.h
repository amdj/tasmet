// system.h, created March 19th 2014
// Author: J.A. de Jong

// A system contains one or more (un)connected (thermoacoustic)
// segments. Some segments require boundary conditions. A System()
// instance combines segments and boundary conditions, such that a
// full error vector can be created, a right hand side and the
// Jacobian.

#pragma once
#ifndef _SYSTEM_H_
#define _SYSTEM_H_
#define MAXNDOFS 20000
#include <memory>
#include <vtypes.h>
#include "tube.h"
#include "bcvertex.h"
#include "globalconf.h"
#include "systemhelpers.h"

#define MAXSEGS 30

namespace tasystem{
  SPOILNAMESPACE

  using segment::Seg;
  using segment::BcVertex;
  using tube::Tube;
  using segment::SegCoupling;

  void coupleSegs(Seg& seg1,Seg& seg2,SegCoupling); // Couple two segments  
  
  
  class TAsystem;
  using segment::Seg;

  class TAsystem{
  private:
    vector<Seg*> segs;
    vector<BcVertex*> bcvertices;
    Globalconf gc;
    // vector<us> startdof;	// Vector containing the starting degree of freedom for segment number # 
    // vector<us> enddof;		// Vector containing the last dof belonging to segment number #
  private:
    us Nsegs=0;			// Number of segments
    us Ndofs=0;
    us Nbc=0;
    bool isInit=false;
    
    
  public:
    TAsystem(const Globalconf& g);
    ~TAsystem();
    TAsystem(const TAsystem& o);
    TAsystem& operator=(const TAsystem& other);
    void Init();
    void show();
    us getNsegs() const {return Nsegs;}
    us getNbc() const {return Nbc;}
    // System with a
    // vector of segments
    vd Error();			// Total error vector
    vd GetRes();			// Extract result vector
    void SetRes(vd resvec);	// Set result vector
    void addseg(const Seg& s);	// Add a segment to the system. It creates a copy.
    void addbc(const BcVertex& vertex);
    BcVertex* getBc(us nr) const;
    // void delseg(us n); // Not yet implementen. Delete a segment from the system (we have to determine how elaborated the API has to be.)
    void setGc(const Globalconf& gc); // Reset globalconf configuration
    dmat Jac();		// Return Jacobian matrix    
    Seg* operator[](us i) const;    
    Seg* getSeg(us i) const {return (*this)[i];} // Easier for cython wrapping
    arma::uvec segfirstcol=zeros<arma::uvec>(MAXSEGS);
    arma::uvec segndofs=zeros<arma::uvec>(MAXSEGS);

  private:
    // A vector of boundary conditions is required
    void CheckInit();
    void cleanup();
    void setnodes(us segnr,us nL,us nR);
    bool hasinit=false;
    // friend void copysegs(TAsystem& to,const TAsystem& from);
  };				// class System
  
} // namespace tasystem

#endif /* _SYSTEM_H_ */
