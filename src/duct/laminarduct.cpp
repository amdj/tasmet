// laminarduct.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////


#include "laminarduct.h"
#include "bccell.h"
#include "soltw.h"

#include "globalconf.h"
#include "geom.h"
#include "var.h"

#include "laminardrag.h"
#include "hopkinsheat.h"

// Tried to keep the method definition a bit in order in which a
// duct is created, including all its components. First a duct is
// created, which has a geometry and a pointer to the global
// configuration. Moreover, the duct has gridpoints,
// "LaminarDuctCell" instants. Of these, a duct has gp of them,
// stored in a vector. In each gridpoint, variables live, which
// represent the current solution. Moreover, we have equations in
// each gridpoint. More precisely, in the final solution the
// continuity, momentum, energy and a suitable equation of state
// should hold.


namespace duct {
  using tasystem::Globalconf;
  using tasystem::TaSystem;
  using tasystem::var;

  LaminarDuct::LaminarDuct(const Geom& geom,d Twl,d Twr):
    LaminarDuct(geom)
  {
    // DO SOMETHING with Tinit!
    if(Twr<0)
      Twr=Twl;
    d L=geom.L();
    const vd& vx=geom.vx_vec();
    Tinit=Twl+(Twr-Twl)*vx/L;
  }
  LaminarDuct::LaminarDuct(const LaminarDuct& o,const TaSystem& sys):
    Duct(o,sys),
    insulated(o.insulated),
    Tinit(o.Tinit)
  {
    laminardrag=new drag::LaminarDragResistance(*this);
    hopkinsheat=new HopkinsHeatSource(*this);
    TRACE(15,"LaminarDuct::LaminarDuct(copy)");
    if(Tinit.size()>0){
      assert(cells.size()>0);
      for(auto& cell : cells){
	tasystem::var Tvar(*gc,Tinit(cell->i));
	cell->setResVar(Varnr::Ts,Tvar);      
	cell->setResVar(Varnr::Tw,Tvar);      
	cell->setResVar(Varnr::T,Tvar);
      }  // // Set time-avg data to make solving bit easier
    }
  }
  void LaminarDuct::init() {
    TRACE(45,"LaminarDuct::init()");
    Duct::init();
  }
  void LaminarDuct::show(us s) const {
    TRACE(15,"LaminarDuct::show()");
    cout << "Laminar Duct\n";
    Duct::show(s);
  }
  LaminarDuct::LaminarDuct(const Geom& geom):
    Duct(geom)
  {
    // Fill vector of gridpoints with data:
    TRACE(13,"LaminarDuct constructor()...");
  }
  void LaminarDuct::setInsulated(bool ins){
    TRACE(15,"LaminarDuct::setInsulated()");
    if(ins && hasSolid())
      throw MyError("When a solid is present such as in a Stack type of duct,"
		    " the Duct cannot be set to insulate.");
    insulated=ins;
  }
  void LaminarDuct::setVarsEqs(Cell& c) const {
    TRACE(15,"LaminarDuct::setVarsEqs()");
    Duct::setVarsEqs(c);

    // Both when insulated as when a solid is present (in case *this
    // is of type Stack), the wall temperature is a dependent
    // variable. The SolTw equation deals with what equation it need
    // to solve in which case.
    if(isInsulated() || hasSolid()){

      auto& eqs=c.getEqs();
      auto& vars=c.getVars();

      // LaminarDuct is not a friend of Cell. Therefore, a const_cast
      // hack is used to put these variables in the vector of
      // variables to solve for.
      vars.push_back(&const_cast<var&>(c.Tw()));
      eqs.insert({EqType::SolTwEq,new SolTw(c,*this)});
    }	// isInsulated() || hasSolid()
  }
  
  LaminarDuct::~LaminarDuct(){
    delete laminardrag;
    delete hopkinsheat;
  }
  
} /* namespace duct */

//////////////////////////////////////////////////////////////////////
