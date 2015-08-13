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


namespace tube {

  #ifdef SWIG
  %catches(std::exception,...) ConnectorVolume::ConnectorVolume(d volume);
  #endif // SWIG
  
  #ifndef SWIG
  struct Connection{
    us segnr;
    Pos position;
    Connection(us segnr,Pos position):segnr(segnr),position(position){}
  };
  class Tube;
  class BcCell;
  struct TubeConnection:public Connection{
    const Tube* t=nullptr;
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

    // Vector containing all Tube connections
    std::vector<TubeConnection> tubeConnections;
    std::vector<PistonConnection> pistonConnections;
    
    tasystem::var p_,T_,rho_;
    ConnectorVolume(const tasystem::TaSystem&,const ConnectorVolume& other);
  public:
    ConnectorVolume(const ConnectorVolume& other)=delete;
    ConnectorVolume(d volume);
    ~ConnectorVolume();
    segment::Seg* copy(const tasystem::TaSystem& s) const {
      return new ConnectorVolume(s,*this);
    }
    const tasystem::var& p() const {return p_;}
    const tasystem::var& rho() const {return rho_;}
    const tasystem::var& T() const {return T_;}

    void addTube(us segnr,Pos position);
    void addPiston(us segnr,Pos position);

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
    virtual vd getRes() const; // Get a result vector
    virtual void domg(vd&) const;	// Derivative of error w.r.t. base frequency.
    virtual void setRes(const vd& res);  // Setting result vector
    virtual void dmtotdx(vd&) const; // Derivative of current fluid mass in
    // system to all dofs.
    #endif

  };
  
} // namespace tube

#endif // CONNECTORVOLUME_H
//////////////////////////////////////////////////////////////////////
