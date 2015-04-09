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
#include "arma_eigen.h"

#ifndef SWIG
namespace segment{
  class Seg;
  class Connector;
}

namespace tube{
  class Tube;
}
#endif

namespace tasystem{
  #ifndef SWIG
  SPOILNAMESPACE

  typedef Eigen::VectorXd evd;
  typedef Eigen::SparseMatrix<double> esdmat;

  class TripletList;
  #endif  

  #ifdef SWIG
  %catches(std::exception,...) TaSystem::getTube(us i) const;
  %catches(std::exception,...) TaSystem::TaSystem();
  %catches(std::exception,...) TaSystem::TaSystem(const Globalconf&);
  %catches(std::exception,...) TaSystem::Error();
  %catches(std::exception,...) TaSystem::init();
  %catches(std::exception,...) TaSystem::showJac();
  %catches(std::exception,...) TaSystem::show(us detailnr=0);
  %catches(std::exception,...) TaSystem::operator+=(const segment::Seg&);
  %catches(std::exception,...) TaSystem::operator+=(const segment::Connector&);
  #endif // SWIG

  class TaSystem{
  protected:
    bool hasInit=false;
    vector<segment::Seg*> segs;		
    vector<segment::Connector*> connectors;    // Yes, connectors are just like segments
    bool driven=true;
  protected:
    Globalconf gc_;             // Global configuration parameters
  public:
    void setGc(const Globalconf& gc); // Reset globalconf configuration
    const Globalconf& gc() const {return gc_;} // Reset globalconf configuration
    TaSystem():gc_(Globalconf::airSTP(0,100)){}
    TaSystem(const Globalconf& g);
    TaSystem(const TaSystem& o);
    TaSystem& operator=(const TaSystem& other)=delete;

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
    vd Error() {return math_common::EigenToArma(error());}// Total error vector
    virtual void show(us detailnr=0);
    virtual void init();

    #ifndef SWIG
    virtual evd error();			// Total error vector
    virtual evd getRes();			// Extract result vector
    virtual esdmat jac(d dummy=-1);		// Return Jacobian matrix

    void setRes(const evd& res);
    #endif
    virtual void setRes(const vd& resvec);	// Set result vector
    vd GetRes();			// Extract result vector
    void setNf(us);

    // Reset amplitude data in higher harmonics
    void resetHarmonics();
    // void delseg(us n); // Not yet implemented.  Delete a segment
    // from the system (we have to determine how elaborated the API
    // has to be.)
    const tube::Tube& getTube(us i) const;
    us getNDofs() const;	// Compute DOFS in system, set     
    us getNEqs() const;    

    #ifndef SWIG                // The unsafe access methods
    segment::Seg* operator[](us i) const;    
    segment::Seg* getSeg(us i) const; // Easier for cython wrapping
    #endif
    d getCurrentMass();	// Return current mass in system [kg]
    void checkInit();
  protected:
    void jacTriplets(TripletList&);
    void setDofEqNrs();
    // A vector of boundary conditions is required
    void cleanup();

  };				// class System
  
} // namespace tasystem

#endif /* _SYSTEM_H_ */
