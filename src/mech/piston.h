// piston.h
//
// Author: J.A. de Jong 
//
// Description:
// This segment comprises a Piston with two volumes on both sides.
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef PISTON_H
#define PISTON_H
#include "seg.h"
#include "var.h"
#include "constants.h"

// A mechanical-fluid interaction segment
namespace mech {
  
  #ifdef SWIG
  %catches(std::exception,...) PistonConfiguration::PistonConfiguration(d Sl,d Sr,d V0l,d V0r,d M,d Km,d Cm,d Stl=-1,d Str=-1);
  %catches(std::exception,...) Piston::Piston(const PistonConfiguration& pc,bool arbitrateMass=false);
  #endif // SWIG

  struct PistonConfiguration{
    PistonConfiguration(d Sl,d Sr,d V0l,d V0r,d M,d Km,d Cm,d Stl=-1,d Str=-1);
    ~PistonConfiguration(){}
    d M;                        // The piston mass [kg]
    d Sr;                        // The piston right area [m^2]
    d Sl;                        // The left area [m^2]

    d Km;                        // The mechanical spring stiffness
    // [N/m]
    d Cm;                       // Mechanical damping constant

    d V0l,V0r;                   // Front and back volume [m^3]

    // Total cross-sectional area of fluid in contact with wall.
    d Stl=-1,Str=-1;
    
  };

  class Piston: public segment::Seg {
    us firsteqnr;               // First equation of this segment
    PistonConfiguration pc;
    // These booleans determine what kind of equations to solve.
    mutable bool leftConnected=false,rightConnected=false;
    d massL=-1,massR=-1;

    tasystem::var xp_; // Piston position, front
    tasystem::var Fp_;                   // Force on piston (externally
    // applied)
    // ml: mass flow out of left volume
    // mr: mass flow out of right volume
    // 
    tasystem::var pl_,pr_,rhol_,rhor_,Tl_,Tr_,ml_,mr_,mHl_,mHr_;

    // Prescribed mean temperature in the volumes. Can be set using
    // setT0().
    d T0=-1;
    bool arbitrateMass=false;
    
    #ifndef SWIG
    Piston(const Piston& other)=delete;
    Piston& operator=(const Piston& other)=delete;    
    #endif // ifndef SWIG
    Piston(const tasystem::TaSystem&,const Piston& other);
  public:
    Piston(const PistonConfiguration& pc,bool arbitrateMass=false);
    ~Piston();
    segment::Seg* copy(const tasystem::TaSystem& s) const {
      return new Piston(s,*this);
    }
    void setT0(d T01){T0=T01;}
    d getT0() const {return T0;}
    const PistonConfiguration& getPc() const {return pc;}
    const tasystem::var& Fpiston() const { return Fp_;}
    const tasystem::var& xpiston() const { return xp_;}
    tasystem::var upiston() const;
    const tasystem::var& pl() const {return pl_;}
    const tasystem::var& pr() const {return pr_;}
    const tasystem::var& p(Pos p) const {return p==Pos::right?pr_:pl_;}
    const tasystem::var& m(Pos p) const {return p==Pos::right?mr_:ml_;}
    const tasystem::var& T(Pos p) const {return p==Pos::right?Tr_:Tl_;}
    const tasystem::var& mH(Pos p) const {return p==Pos::right?mHr_:mHl_;}
    const tasystem::var& rhol() const {return rhol_;}
    const tasystem::var& rhor() const {return rhor_;}
    const tasystem::var& Tl() const {return Tl_;}
    const tasystem::var& Tr() const {return Tr_;}
    int arbitrateMassEq() const;
    // Post-processing the left and right volume (vars)
    tasystem::var Vl() const;
    tasystem::var Vr() const;

    virtual vd error() const;

    #ifndef SWIG
    // If a side of the piston is not connected, different equations
    // are solved in that part. See documentation for details
    void setConnected(Pos pos) const;
    bool isConnected(Pos pos) const;

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
  inline const Piston& asPiston(const segment::Seg& s){
    return dynamic_cast<const Piston&>(s);
  }

  
} // namespace mech

#endif // PISTON_H
//////////////////////////////////////////////////////////////////////
