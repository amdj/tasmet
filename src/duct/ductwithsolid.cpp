// ductwithsolid.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////

#include "ductwithsolid.h"
#include "solidenergy.h"
#include "solid.h"
#include "bccell.h"
#include "var.h"
#include "geom.h"

namespace duct {
  using tasystem::var;

  DuctWithSolid::DuctWithSolid(const string& solidstr,d ksfrac):
    ksfrac(ksfrac)
  {
    TRACE(15,"DuctWithSolid::DuctWithSolid()");
    solid = new solids::Solid(solidstr);    
  }
  DuctWithSolid::DuctWithSolid(const DuctWithSolid& o):
    Qsin(o.Qsin),
    ksfrac(o.ksfrac)
  {
    TRACE(15,"DuctWithSolid::DuctWithSolid()");
    solid = new solids::Solid(*o.solid);
  }
  void DuctWithSolid::setVarsEqs(Cell& c) const {
    TRACE(15,"DuctWithSolid::setVarsEqs()");
    auto& eqs=c.getEqs();
    auto& vars=c.getVars();

    // Add solid temperature as variable
    vars.push_back(&const_cast<var&>(c.Ts()));

    // Add solid energy equation as equation
    eqs.insert({EqType::Sol,new SolidEnergy(c,solid,ksfrac)});

    d Lduct=geom().L();
    d Lcell=c.xr-c.xl;
    d Qsincell=Qsin*Lcell/Lduct;
    static_cast<SolidEnergy*>(c.getEqs().at(EqType::Sol))->setQin(Qsincell);

    if((!c.left()) || (!c.right())) {
      // Downcast as we know the cells at the boundaries are BcCells
      auto& d=static_cast<BcCell&>(c);
      vars.push_back(&const_cast<var&>(d.Tsbc()));
    }

  }
  DuctWithSolid::~DuctWithSolid(){
    delete solid;
  }
  void DuctWithSolid::setQsin(d Qsin) {
    TRACE(15,"DuctWithSolid::setQsin()");
    this->Qsin=Qsin;
  }

} // namespace duct
//////////////////////////////////////////////////////////////////////
