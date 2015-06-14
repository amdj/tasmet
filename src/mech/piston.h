// piston.h
//
// Author: J.A. de Jong 
//
// Description:
//
//////////////////////////////////////////////////////////////////////
#pragma once
#ifndef PISTON_H
#define PISTON_H
#include "seg.h"
#include "var.h"
#include "constants.h"

// A mechanical-fluid interaction segment
namespace mech {
  
  #ifndef SWIG
  class PistonVolumeEq;
  #endif

  class Piston: public segment::Seg {
    us firsteqnr;               // First equation of this segment
    // These booleans determine what kind of equations to solve.
    bool leftConnected=false,rightConnected=false;
    d M;                        // The piston mass [kg]
    d Sr;                        // The piston right area [m^2]
    d Sl;                        // The left area [m^2]

    d Stl,Str;                  // Total cross-sectional area of fluid
                                // in contact with wall.
    d Km;                        // The mechanical spring stiffness
    // [N/m]
    d Cm;                       // Mechanical damping constant

    d V0l,V0r;                   // Front and back volume [m^3]
    d massL=-1,massR=-1;

    tasystem::var xp_; // Piston position, front
    tasystem::var Fp_;                   // Force on piston (externally
                                        // applied)

    // Left pressure,right pressure

    // ml: mass flow out of left volume
    // mr: mass flow out of right volume
    tasystem::var pl_,pr_,rhol_,rhor_,Tl_,Tr_,ml_,mr_;
    
    Piston(const tasystem::TaSystem&,const Piston& other);
  public:
    Piston(const Piston& other)=delete;
    Piston(d Sl,d Sr,d V0l,d V0r,d M,d Km,d Cm,d Str=-1,d Srl=-1,
           bool leftConnected=false,bool rightConnected=false):
      Seg(),
      leftConnected(leftConnected),rightConnected(rightConnected),
      M(M),Sr(Sr),Sl(Sl),Km(Km),Cm(Cm),V0l(V0l),V0r(V0r)
    {
      TRACE(15,"Piston()");
    }
    ~Piston();
    segment::Seg* copy(const tasystem::TaSystem& s) const {
      return new Piston(s,*this);
    }

    const tasystem::var& Fpiston() const { return Fp_;}
    const tasystem::var& xpiston() const { return xp_;}
    const tasystem::var& pl() const {return pl_;}
    const tasystem::var& pr() const {return pr_;}
    const tasystem::var& rhol() const {return rhol_;}
    const tasystem::var& rhor() const {return rhor_;}
    const tasystem::var& Tl() const {return Tl_;}
    const tasystem::var& Tr() const {return Tr_;}
    #ifndef SWIG
    // If a side of the piston is not connected, different equations
    // are solved in that part. See documentation for details
    void setConnected(Pos pos,bool con);
    bool isConnected(Pos pos,bool con) const;

    virtual vd error() const;
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
  
} // namespace mech

#endif // PISTON_H
//////////////////////////////////////////////////////////////////////
