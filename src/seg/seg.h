#pragma once
#ifndef _SEG_H_
#define _SEG_H_

#include <memory>
#include "vtypes.h"
#include "segconbase.h"


namespace segment{

  class Seg:public SegConBase{
  public:
    Seg();
    Seg(const Seg& o);
    virtual ~Seg(){}            // We do not own the gc instance
    bool operator==(const Seg& other) const {return (this==&other);}
    // Pure virtual functions
    const tasystem::Globalconf& getGc() const;
    virtual Seg* copy() const=0;

    virtual void resetHarmonics()=0;
    virtual std::string getType() const=0;		   // This param is
    // important for connecting the segments
    // ------------------------------ config methods
    virtual void setDofNrs(us firstdofnr)=0;
    virtual us getNDofs() const=0;
    virtual d getCurrentMass() const=0;
    // ------------------------------
    virtual void domg(vd&) const=0;	// Derivative of error w.r.t. base frequency.
    virtual void setRes(const vd& res)=0;  // Setting result vector
    virtual void setRes(const Seg& res)=0; // Copying contents
    virtual void dmtotdx(vd&) const=0; // Derivative of current fluid mass in
    // system to all dofs.
    virtual vd getRes() const=0; // Get a result vector
  };
  
} // Namespace segment


#endif /* _SEG_H_ */



