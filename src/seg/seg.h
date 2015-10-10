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
    #ifndef SWIG
    Seg(const Seg& o)=delete;
    Seg& operator=(const Seg&)=delete;
    #endif // ifndef SWIG
  public:
    virtual segment::Seg* copy(const tasystem::TaSystem&) const=0;
    virtual ~Seg(){}            // We do not own the gc instance

    // This function can be used to set a phase constraint for this
    // segment. This default implementation throws an error telling
    // that the segment is unable to provide a DOF to constrain the
    // phase on
    virtual void setPhaseContraint(tasystem::PhaseConstraint);
    #ifndef SWIG
    // Pure virtual functions

    // Dof number for which an additional constraint can be set to
    // uniquely determine the phase of the system. (For EngineSystems,
    // the phase of higher harmonics is arbitrary and needs to be
    // fixed). This method returns the Dof on which an EngineSystem
    // puts a phase constraint as extra equation when the extra DOF
    // omega is introduced.
    virtual int providePhaseDof() const {return -1;}
    virtual d phaseDofValue() const {return 0;}
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



