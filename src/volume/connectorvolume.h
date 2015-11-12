// connectorvolume.h
//
// Author: J.A. de Jong 
//
// Description:
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef CONNECTORVOLUME_H
#define CONNECTORVOLUME_H

#include "seg.h"
#include "var.h"
#include "constants.h"

#ifdef SWIG
%catches(std::exception,...) duct::asConnnectorVolume(const segment::Seg&);
#endif // SWIG

namespace duct {

  #ifdef SWIG
  %catches(std::exception,...) ConnectorVolume::ConnectorVolume(d volume);
  #endif // SWIG
  
  #ifndef SWIG
  // Helper struct which is used to contain the ducts and pistons
  // which are connected to a ConnectorVolume
  struct Connection{
    std::string segid;
    Pos position;
    Connection(const string& segid,Pos position):segid(segid),position(position){}
  };
  class Duct;
  class BcCell;

  struct DuctConnection:public Connection{
    const Duct* t=nullptr;
    const BcCell* c=nullptr;
    // As we only need the base constructor
    using Connection::Connection;
    void setPtr(const tasystem::TaSystem&);
  };

  struct PistonConnection:public Connection{};
  #endif

  class ConnectorVolume: public segment::Seg {
    us firsteqnr;               // First equation of this segment
    d volume;                   // The fluid volume
    // Vector containing all Duct connections
    std::vector<DuctConnection> ductConnections;
    std::vector<PistonConnection> pistonConnections;

    bool arbitrateMass=false;	// Whether this segment arbitrates the
				// mass in the total system.
    
    tasystem::var p_,T_,rho_;
    ConnectorVolume(const tasystem::TaSystem&,const ConnectorVolume& other);
    #ifndef SWIG
    ConnectorVolume(const ConnectorVolume& other)=delete;
    ConnectorVolume& operator=(const ConnectorVolume& other)=delete;    
    #endif // ifndef SWIG


  public:
    ConnectorVolume(d volume,bool arbitrateMass=false);
    ~ConnectorVolume();
    segment::Seg* copy(const tasystem::TaSystem& s) const {
      return new ConnectorVolume(s,*this);
    }
    int arbitrateMassEq() const;
    const tasystem::var& p() const {return p_;}
    const tasystem::var& rho() const {return rho_;}
    const tasystem::var& T() const {return T_;}

    void addDuct(const string& segid,Pos position);
    void addPiston(const string& segid,Pos position);

    virtual vd error() const;

    #ifndef SWIG

    virtual void setEqNrs(us firstdofnr);    
    virtual us getNEqs() const;
    virtual void show(us) const;
    virtual void jac(tasystem::Jacobian&) const;
    virtual void updateNf();  // Update nr of frequencies.

    virtual void resetHarmonics();
    virtual void setDofNrs(us firstdofnr);
    virtual us getNDofs() const;
    virtual d getMass() const;
    // ------------------------------
    #endif
    virtual vd getRes() const; // Get a result vector
    #ifndef SWIG
    virtual void domg(vd&) const;	// Derivative of error w.r.t. base frequency.
    virtual void setRes(const vd& res);  // Setting result vector
    virtual void dmtotdx(vd&) const; // Derivative of current fluid mass in
    // system to all dofs.
    #endif

  };
  
  inline const ConnectorVolume& asConnectorVolume(const segment::Seg& s){
    return dynamic_cast<const ConnectorVolume&>(s);
  }

} // namespace duct

#endif // CONNECTORVOLUME_H
//////////////////////////////////////////////////////////////////////
