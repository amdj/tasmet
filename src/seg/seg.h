#pragma once
#ifndef _SEG_H_
#define _SEG_H_

#include <memory>
#include "vtypes.h"
#include "segconbase.h"

namespace segment{

  class Seg:public SegConBase{
  protected:
    Seg();
    Seg(const Seg& o,const tasystem::TaSystem& s): SegConBase(o,s){}
    Seg(const Seg& o)=delete;
    Seg& operator=(const Seg&)=delete;
  public:
    virtual segment::Seg* copy(const tasystem::TaSystem&) const=0;
    virtual ~Seg(){}            // We do not own the gc instance
    #ifndef SWIG
    // Pure virtual functions
    virtual void resetHarmonics()=0;
    // important for connecting the segments
    // ------------------------------ config methods
    virtual void setDofNrs(us firstdofnr)=0;
    virtual us getNDofs() const=0;
    virtual d getMass() const=0;
    // ------------------------------
    virtual vd getRes() const=0; // Get a result vector
    virtual void domg(vd&) const=0;	// Derivative of error w.r.t. base frequency.
    virtual void setRes(const vd& res)=0;  // Setting result vector
    virtual void dmtotdx(vd&) const=0; // Derivative of current fluid mass in
    // system to all dofs.
    #endif

  };

} // Namespace segment


#endif /* _SEG_H_ */



