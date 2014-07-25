#pragma once
#ifndef _SYSTEMHELPERS_H_
#define _SYSTEMHELPERS_H_
#include "seg.h"
#include "bcvertex.h"




namespace tasystem{

  using segment::Seg;
  using segment::Vertex;
  using segment::BcVertex;
  using segment::connectpos;

  class TAsystem;
  
  void copyallsegsbc(TAsystem& to,const TAsystem& from);  
  void connectbc(Seg&,const BcVertex&);
  Seg* copyseg(const Seg& orig);
  BcVertex* copybc(const BcVertex& orig);
  Vertex* vertexfrombc(BcVertex*);

} // namespace tasystem

#endif /* _SYSTEMHELPERS_H_ */




