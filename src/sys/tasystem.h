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

#include <memory>
#include "vtypes.h"
#include "globalconf.h"

#ifndef SWIG
namespace segment{
  class Seg;
  class Connector;
}
namespace mech {
  class Piston;
} // namespace mech

namespace tube{
  class Tube;
}
#endif

namespace tasystem{
  #ifndef SWIG
  SPOILNAMESPACE

  class TripletList;
  #endif  

  #ifdef SWIG
  %catches(std::exception,...) TaSystem::getTube(us i) const;
  %catches(std::exception,...) TaSystem::getPiston(us i) const;
  %catches(std::exception,...) TaSystem::getSeg(us i) const;
  %catches(std::exception,...) TaSystem::getConnector(us i) const;
  %catches(std::exception,...) TaSystem::TaSystem();
  %catches(std::exception,...) TaSystem::TaSystem(const Globalconf&);
  %catches(std::exception,...) TaSystem::Error();
  %catches(std::exception,...) TaSystem::init();
  %catches(std::exception,...) TaSystem::setRes(const vd& res);
  %catches(std::exception,...) TaSystem::showJac();
  %catches(std::exception,...) TaSystem::show(us detailnr=0);
  %catches(std::exception,...) TaSystem::operator+=(const segment::Connector&);
  %catches(std::exception,...) TaSystem::operator+=(const segment::Seg&);
  #endif // SWIG

  class TaSystem{
  protected:
    bool hasInit=false;
    // This is the mass in the sytem. 
    d mass_=-1;
    // Tells the TaSystem which Dof should be overwritten with the
    // mass arbitration equation. The special value of -1 means, that
    // mass is not arbitrated. This can be the case if for example a
    // pressure boundary condition is present.
    int arbitrateMassEq=-1;

    vector<segment::Seg*> segs;		
    vector<segment::Connector*> connectors;    // Yes, connectors are just like segments
    bool driven=true;
    Globalconf gc_;             // Global configuration parameters
  public:
    const Globalconf& gc() const {return gc_;} // Reset globalconf configuration
    TaSystem():gc_(Globalconf::airSTP(0,100)){}
    TaSystem(const Globalconf& g);
    TaSystem(const TaSystem& o);
    TaSystem& operator=(const TaSystem& other)=delete;


    // Set globalconf configuration. Applies updateNf as well.
    void setGc(const Globalconf& gc);
    

    // Set and get the mass in the system. If the mass is not set
    // before initializing, the mass is computed from the segment's
    // intial configuration.
    void setMass(d mass){mass_=mass;}
    d getMass() const {return mass_;}

    virtual ~TaSystem();
    virtual TaSystem* copy() const {return new TaSystem(*this);}

    void setDriven(bool dr) {driven=dr;}
    bool isDriven() const {return driven;}
    us nSegs() const {return segs.size();}
    us nConnectors() const {return connectors.size();}
    TaSystem& operator+=(const segment::Connector& c);
    TaSystem& operator+=(const segment::Seg& s);	// Add a segment to the
    // system. It creates a copy

    dmat showJac();
    virtual void show(us detailnr=0);
    virtual void init();
    void checkInit();
    virtual vd Error();			// Total error vector
    virtual vd getRes();			// Extract result vector
    virtual void setRes(const vd& resvec);	// Set result vector
    #ifndef SWIG
    // Compute Jacobian matrix. The dampfac value is used in an
    // EngineSystem
    arma::sp_mat jac(d dampfac=1);		// Return Jacobian matrix
    #endif

    // Change Nf in the system, while keeping the results.
    void updateNf(us);
    // Reset amplitude data in higher harmonics
    void resetHarmonics();
    // void delseg(us n); // Not yet implemented.  Delete a segment
    // from the system (we have to determine how elaborated the API
    // has to be.)
    const tube::Tube& getTube(us i) const;
    const mech::Piston& getPiston(us i) const;
    us getNDofs() const;	// Compute DOFS in system, set     
    us getNEqs() const;    

    const segment::Connector* getConnector(us i) const;    
    const segment::Seg* getSeg(us i) const;

  protected:
    // Return the mass of the system. If mass is arbitrated, this
    // should be equal to the result of getMass if the solution is
    // converged.
    d getCurrentMass();
    vd dmtotdx() const;         // Derivative of total mass to DOF x
    virtual TripletList jacTriplets(d dummy);
    virtual void cleanup();
  };				// class System
  
} // namespace tasystem

#endif /* _SYSTEM_H_ */

