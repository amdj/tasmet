// laminarduct.cpp
//
// last-edit-by: J.A. de Jong 
// 
// Description:
//
//////////////////////////////////////////////////////////////////////


#include "laminarduct.h"
#include "bccell.h"
#include "solid.h"
#include "solidenergy.h"
#include "soltw.h"


#include "globalconf.h"
#include "geom.h"
#include "var.h"
// Tried to keep the method definition a bit in order in which a
// tube is created, including all its components. First a tube is
// created, which has a geometry and a pointer to the global
// configuration. Moreover, the tube has gridpoints,
// "LaminarDuctCell" instants. Of these, a tube has gp of them,
// stored in a vector. In each gridpoint, variables live, which
// represent the current solution. Moreover, we have equations in
// each gridpoint. More precisely, in the final solution the
// continuity, momentum, energy and a suitable equation of state
// should hold.
using tasystem::Globalconf;
using tasystem::TaSystem;
using tasystem::var;

namespace tube {
  LaminarDuct::LaminarDuct(const Geom& geom,const string& solid):
    LaminarDuct(geom)
  {
    setSolid(solid);
  }

  LaminarDuct::LaminarDuct(const Geom& geom,d Tl,d Tr):
    LaminarDuct(geom)
  {
    this->Tl=Tl;
    this->Tr=Tr;
  }
  LaminarDuct::LaminarDuct(const LaminarDuct& o,const TaSystem& sys):
    Tube(o,sys),
    laminardrag(*this),
    hopkinsheat(*this),
    Tl(o.Tl),
    Tr(o.Tr)
  {
    TRACE(15,"LaminarDuct::LaminarDuct(copy)");
    if(o.solid)
      solid=new solids::Solid(*o.solid);
  }
  void LaminarDuct::setSolid(const string& solidname,d ksfrac) {
    TRACE(15,"LaminarDuct::setSolid()");
    if(solid)			// Delete old solid
      delete solid;
    solid=new solids::Solid(solidname);
  }
  void LaminarDuct::init() {
    TRACE(45,"LaminarDuct::init()");
    Tube::init();
    d L=geom().L();
    if(Tl<0)
      Tl=gc->T0();
    if(Tr<0)
      Tr=gc->T0();
    d T;
    assert(cells.size()>0);
    const vd& vx=geom().vx_vec();
    // Tmirror=Tl+(Tr-Tl)*math_common::skewsine(xv/L);
    vd Tmirror=Tl+(Tr-Tl)*vx/L;
    for(auto& cell : cells){
      tasystem::var Tvar(*gc,Tmirror(cell->i));
      cell->setResVar(Varnr::Ts,Tvar);      
      cell->setResVar(Varnr::Tw,Tvar);      
      cell->setResVar(Varnr::T,Tvar);
    }  // // Set time-avg data to make solving bit easier

  }
  void LaminarDuct::show(us s) const {
    TRACE(15,"LaminarDuct::show()");
    Tube::show(s);
    string sol;
    if(solid)
      sol=*solid;		// Operator string converts solid to a
				// certain type name.
    else
      sol="none";
    cout << "Solid in system: " << sol << "\n";
  }
  LaminarDuct::LaminarDuct(const Geom& geom):
    Tube(geom),
    laminardrag(*this),
    hopkinsheat(*this)
  {
    // Fill vector of gridpoints with data:
    TRACE(13,"LaminarDuct constructor()...");
  }
  const solids::Solid& LaminarDuct::getSolid() const{
    if(hasSolid())
      return *solid;
    else
      throw MyError("LaminarDuct does not have a solid!");
  }
  void LaminarDuct::setVarsEqs(Cell& c) const {
    TRACE(15,"LaminarDuct::setVarsEqs()");
    Tube::setVarsEqs(c);
    if(hasSolid()) {
      auto& eqs=c.getEqs();
      auto& vars=c.getVars();

      // LaminarDuct is not a friend of Cell. Therefore, a const_cast
      // hack is used to put these variables in the vector of
      // variables to solve for.

      vars.push_back(&const_cast<var&>(c.Ts()));
      vars.push_back(&const_cast<var&>(c.Tw()));
      eqs.insert({EqType::Sol,new SolidEnergy(c,solid)});
      eqs.insert({EqType::SolTwEq,new SolTw(c,*this)});

      if((!c.left()) || (!c.right())) {
	// Downcast as we know the cells at the boundaries are BcCells
	auto& d=static_cast<BcCell&>(c);
	vars.push_back(&const_cast<var&>(d.Tsbc()));
      }
    } // hasSolid()
  }
  LaminarDuct::~LaminarDuct(){
    TRACE(15,"~LaminarDuct()");
    delete solid;
  }
  // vd LaminarDuct::dragCoefVec(us freqnr) const{
  //   TRACE(15,"LaminarDuct::drag_vec()");
  //   vd dragcoef(getNCells());
  //   var drag_varcoef(*gc);
  //   for(us i=0;i<dragcoef.size();i++){
  //     const Cell& cell=getCell(i);
  //     drag_varcoef.setadata(laminardrag.ComplexResistancecoef(cell));
  //     dragcoef(i)=drag_varcoef(freqnr);
  //   }
  //   return dragcoef;
  // }
  
} /* namespace tube */

//////////////////////////////////////////////////////////////////////
