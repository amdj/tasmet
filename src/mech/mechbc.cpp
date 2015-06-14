// mechbc.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////

#include "mechbc.h"
#include "prescribeqty.h"
#include "exception.h"
#include "piston.h"
#include "tasystem.h"
#include "jacobian.h"

#define fDFT (gc->fDFT)
#define iDFT (gc->iDFT)
#define DDTfd (gc->DDTfd)

#define Ns (gc->Ns())
#define Nf (gc->Nf())
#define eye (eye(Ns,Ns))

namespace mech {
  SPOILNAMESPACE  
  using tasystem::JacRow;
  using tasystem::JacCol;

  MechBc::MechBc(us segnr,Varnr var,const tasystem::var& bc):
    Connector(),
    segnr(segnr),
    type(var),
    bc(bc)
  {
    TRACE(15,"MechBc::MechBc()");

    if((type != Varnr::F) && (type!= Varnr::x) && (type!= Varnr::Z)) {
      WARN("Given: " << toString(type));
      throw MyError("Boundary condition can only be set to "
                    "Force, piston displacement, or impedance.");
    }
  }
  MechBc::MechBc(const MechBc& other,const tasystem::TaSystem& sys):
    Connector(other,sys),
    segnr(other.segnr),
    type(other.type),
    bc(other.bc)
  {
    TRACE(15,"MechBc::MechBc()");
    try{
      p=&dynamic_cast<const Piston&>(*sys.getSeg(segnr));
    }
    catch(...) {
      throw MyError("Segment is not a Piston or does not exist.");
    }
    
    setInit(true);
  }

  MechBc::~MechBc(){  }
  us MechBc::getNEqs() const{
    TRACE(15,"us MechBc::getNEqs()");
    return Ns;
  }
  void MechBc::show(us detailnr) const{
    TRACE(15,"void MechBc::show()");
    const char* bctext;
    assert((type==Varnr::F) || (type==Varnr::x) || (type==Varnr::Z));
    if(type==Varnr::F)
      bctext="force b.c."; 
    else if(type==Varnr::x)
      bctext="piston displacement b.c.";
    else
      bctext="Mechanical impedance b.c.";
    cout << "Mechanical domain boundary condition set for segment "<< segnr \
         << ".\n Type of boundary condition:" << bctext << endl;
    if(detailnr>2)
      cout << "Value of bc: " << bc() << endl;
  }
  vd MechBc::error() const{
    TRACE(15,"vd MechBc::error()");

    if(type==Varnr::Z){
      return p->Fpiston()()-bc.freqMultiplyMat()*p->xpiston()();
    }
    else if(type==Varnr::x)
      return p->xpiston()()-bc();
    else
      return p->Fpiston()()-bc();
  }
  void MechBc::jac(tasystem::Jacobian& jac) const{
    TRACE(15,"void MechBc::jac()");
    TRACE(15,"forsteqnr: "<< firsteqnr);
    JacRow bcjac(firsteqnr,2);
    if(type==Varnr::Z){
      bcjac+=JacCol(p->Fpiston(),eye);
      bcjac+=JacCol(p->xpiston(),-bc.freqMultiplyMat());
    }
    else if(type==Varnr::x){
      bcjac+=JacCol(p->xpiston(),eye);      
    }
    else {
      bcjac+=JacCol(p->Fpiston(),eye);      
    }
    jac+=bcjac;
  }
  void MechBc::updateNf(){
    TRACE(15,"void MechBc::updateNf()");
    bc.updateNf();
  }
  
} // namespace mech

//////////////////////////////////////////////////////////////////////
