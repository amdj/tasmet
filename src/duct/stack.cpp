// stack.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////
#include "stack.h"
#include "solid.h"

namespace duct {
  using tasystem::Globalconf;
  using tasystem::TaSystem;

  Stack::Stack(const Geom& geom,const string& solidstr):
    // Should be called as Duct is a
    // virtual base class of LaminardDuct
    // and DuctWithSolid
    Duct(geom),
    // 
    DuctWithSolid(geom,solidstr),
    LaminarDuct(geom)
  {
    TRACE(15,"Stack::Stack()");
    // Initialize a new solid
   
  }
  Stack::Stack(const Stack& o,const TaSystem& sys):
    Duct(o,sys),
    DuctWithSolid(o,sys),
    LaminarDuct(o,sys)
  {
    TRACE(15,"Stack::Stack(copy)");

  }
  void Stack::show(us s) const{
    TRACE(15,"Stack::show()");
    LaminarDuct::show(s);
    cout << "Solid in system: " << string(getSolid()) << "\n";
  }
  void Stack::setVarsEqs(Cell& c) const{
    TRACE(15,"Stack::setVarsEqs()");

    LaminarDuct::setVarsEqs(c);
    DuctWithSolid::setVarsEqs(c);
    
  }
  Stack::~Stack(){ }

} // namespace duct

//////////////////////////////////////////////////////////////////////
