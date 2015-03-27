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
#include "segconnection.h"
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

  class TaSystem{
  protected:
    bool hasInit=false;
    vector<segment::Seg*> segs;		
    vector<segment::Connector*> connectors;    // Yes, connectors are just like segments
    bool driven=true;
  public:
    Globalconf gc;    
    void setGc(const Globalconf& gc); // Reset globalconf configuration

    TaSystem():gc(Globalconf::airSTP(0,100)){}
    TaSystem(const Globalconf& g);
    TaSystem(const TaSystem& o);
    #ifndef SWIG
    TaSystem& operator=(const TaSystem& other);
    #endif
    virtual ~TaSystem();
    virtual TaSystem* copy() const {return new TaSystem(*this);}
    virtual bool init();
    void setDriven(bool dr) {driven=dr;}
    bool isDriven() const {return driven;}
    us nSegs() const {return segs.size();}
    us nConnectors() const {return connectors.size();}
    TaSystem& operator+=(const segment::Connector& c);
    TaSystem& operator+=(const segment::Seg& s);	// Add a segment to the
    // system. It creates a copy

    void showJac(bool force=false);
    vd Error() {return math_common::EigenToArma(error());}// Total error vector
    virtual void show(us detailnr=0);
    #ifndef SWIG
    virtual evd error();			// Total error vector
    virtual evd getRes();			// Extract result vector
    virtual esdmat jac(d dummy=-1);		// Return Jacobian matrix

    void setRes(const evd& res);
    #endif
    virtual void setRes(const vd& resvec);	// Set result vector

    void setNf(us);

    // Reset amplitude data in higher harmonics
    void resetHarmonics();
    // void delseg(us n); // Not yet implemented.  Delete a segment
    // from the system (we have to determine how elaborated the API
    // has to be.)
    tube::Tube* getTube(us i) const throw(std::exception);    
    #ifndef SWIG                // The unsafe access methods
    segment::Seg* operator[](us i) const;    
    segment::Seg* getSeg(us i) const; // Easier for cython wrapping
    #endif
    void setRes(const TaSystem& o);
    d getCurrentMass();	// Return current mass in system [kg]
    bool checkInit(){		// Often called simple method: inline
      if(!hasInit){
        return init(); 
      }
      else
        return hasInit;
    }
  protected:
    us getNDofs() const;	// Compute DOFS in system, set     
    us getNEqs() const;    
    void jacTriplets(TripletList&);
    void setDofEqNrs();
    // A vector of boundary conditions is required
    void copyTaSystem(const TaSystem& other);
    void cleanup();

  };				// class System
  
} // namespace tasystem

#endif /* _SYSTEM_H_ */
