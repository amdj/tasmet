/*
 * linduct.cpp

 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#include "duct.h"
#include "leftcell.h"       // which includes ductbccell.h
#include "rightcell.h"       // which includes ductbccell.h
#include "globalconf.h"
#include "geom.h"
#include "exception.h"
#include "utils.h"

#include "continuity.h"
#include "momentum.h"
#include "mu.h"
#include "energy.h"
#include "state.h"
#include "solidenergy.h"
#include "isentropic.h"


// Tried to keep the method definition a bit in order in which a
// duct is created, including all its components. First a duct is
// created, which has a geometry and a global
// configuration. Moreover, the duct has gridpoints, "Cell"
// instants. Of these, a duct has gp of them, stored in a vector. In
// each gridpoint, variables live, which represent the current
// solution. Moreover, we have equations in each gridpoint. More
// precisely, in the final solution the continuity, momentum, energy
// and a suitable equation of state should hold.

namespace duct {
  using tasystem::PhaseConstraint;
  using tasystem::TaSystem;
  using tasystem::Globalconf;
  using tasystem::Jacobian;


  Duct::Duct(const Geom& geom)
    :Seg(),geom_(geom.copy())
  {
    TRACE(13,"Duct constructor()...");
  }
  Duct::Duct(const Duct& other,const TaSystem& sys):
    Seg(other,sys),
    geom_(other.geom().copy()){
    if(other.pc_)
      this->pc_=new PhaseConstraint(*other.pc_);
    TRACE(13,"Filling vertices. Current size:"<<cells.size());
    // Left *probable* boundary condition
    // cells.emplace_back(new LeftCell(0,g));
    // WARN("Lot wrong here");
    cells.emplace_back(new LeftCell(0,*this));
    for(us i=1;i<getNCells()-1;i++)
      cells.emplace_back(new Cell(i,*this));
    cells.emplace_back(new RightCell(getNCells()-1,*this));

  } // Duct::Duct(copy)
  void Duct::init(){
    TRACE(15,"Duct::init()");

    us nCell=cells.size();    
    assert(nCell==getNCells());
    // And initialize again.
    for(us i=0;i<cells.size();i++){
      TRACE(13,"Starting intialization of Cell "<< i);
      Cell* thiscell=cells[i];
      Cell* left=nullptr;
      Cell* right=nullptr;
      if(i<nCell-1)
        right=cells[i+1];
      if(i>0)
        left=cells[i-1];
      TRACE(15,"Initializing cell");
      thiscell->init(left,right);
    } // for

  }
  void Duct::setVarsEqs(Cell& c) const {
    TRACE(15,"Duct::setVarsEqs()");
    
    auto& vars=c.getVars();
    // Assuming this is the first function that fills the vars vector
    // These do not be deleted!
    vars.clear();
    vars.reserve(constants::nvars_reserve);

    // Assuming this is the first function that fills the eqs map
    auto& eqs=c.getEqs();
    // Equations need to be deleted.
    utils::purge(eqs);

    // These should always be put in
    vars.push_back(&c.rho_);
    vars.push_back(&c.ml_);
    vars.push_back(&c.T_);
    vars.push_back(&c.p_);
    vars.push_back(&c.mu_);

    if((!c.left()) || (!c.right())) {
      auto& d=static_cast<BcCell&>(c);
      vars.push_back(&d.rhobc_);
      vars.push_back(&d.Tbc_);
      vars.push_back(&d.pbc_);
      vars.push_back(&d.ubc_);
      vars.push_back(&d.mHbc_);

      eqs.insert({EqType::BcEqP,new ExtrapolatePressure(d)});
      eqs.insert({EqType::BcEqStateBc,new StateBc(d)});    
      eqs.insert({EqType::BcEqu,new BcVelocity(d)});
    }

    eqs.insert({EqType::Con,new Continuity(c)});
    if(c.left())
      eqs.insert({EqType::Mom,new Momentum(c)});

    eqs.insert({EqType::Ene,new Energy(c)});
    eqs.insert({EqType::Sta,new State(c)});
    eqs.insert({EqType::Mu_is_m_u,new MuEq(c)});    

  }
  Duct::~Duct(){
    TRACE(25,"~Duct()");
    delete geom_;
    delete pc_;
    utils::purge(cells);
  }
  const Cell& Duct::getCell(int i) const{
    if(i>=0)
      return *cells.at(i);
    else
      return *cells.at(cells.size()+i);
  }
  const Cell& Duct::operator[](us i) const{
    return *cells[i];
  }
  us Duct::getNCells() const {return geom().nCells();}
  const Geom& Duct::geom() const {return *geom_;}
  void Duct::setDofNrs(us firstdof){
    TRACE(13,"Duct::setDofNrs("<<firstdof<<")");
    assert(cells.size()>0);
    for(Cell* cell: cells){
      cell->setDofNrs(firstdof);
      firstdof+=cell->getNDofs();
    }
  }
  void Duct::setEqNrs(us firsteq){
    TRACE(13,"Duct::SetEqNrs("<<firsteq<<")");
    assert(cells.size()>0);
    for(auto cell=cells.begin();cell!=cells.end();cell++){
      (*cell)->setEqNrs(firsteq);
      firsteq+=(*cell)->getNEqs();
    }
  }
  
  us Duct::getNDofs() const {
    TRACE(10,"Duct::getNDofs()");
    us ndofs=0;
    for(auto v=cells.begin();v!=cells.end();v++)
      ndofs+=(*v)->getNDofs();
    return ndofs;
  }
  us Duct::getNEqs() const {
    TRACE(10,"Duct::getNEqs()");
    us ndofs=0;
    for(auto v=cells.begin();v!=cells.end();v++)
      ndofs+=(*v)->getNEqs();
    return ndofs;
  }  

  d Duct::getMass() const{
    TRACE(8,"Duct::getMass()");
    assert(cells.size()>0);
    d mass=0;
    for(auto cell=cells.begin();cell!=cells.end();cell++){
      mass+=(*cell)->getMass();
    }
    return mass;
  }
  void Duct::updateNf(){
    TRACE(18,"Duct::updateNf()");
    assert(cells.size()>0);
    for(auto v: cells){
      v->updateNf();
    }
  }
  // vc Duct::heatQ() const {
  //   TRACE(15,"Duct::heatQ()");
  //   vc res(getNCells());
  //   for(us i=0; i< getNCells(); i++) {
  //     res(i)=static_cast<const HopkinsHeatSource&>(getHeatSource()).HeatTransferCoefQ(*cells[i])(1);
  //   }
  //   return res;
  // }
  // vc Duct::heatH() const {
  //   TRACE(15,"Duct::heatQ()");
  //   vc res(getNCells());
  //   for(us i=0; i< getNCells(); i++) {
  //     res(i)=static_cast<const HopkinsHeatSource&>(getHeatSource()).HeatTransferCoefH(*cells[i])(1);
  //   }
  //   return res;
  // }
  const BcCell& Duct::bcCell(Pos p) const{
    TRACE(3,"Duct::bcCell()");
    if(p==Pos::left){
      assert(cells[0]);
      return static_cast<const BcCell&>(*cells[0]);
    }
    else{
      return static_cast<const BcCell&>(**(cells.end()-1));
    }
  }
  vd Duct::getx() const {
    TRACE(10,"Duct::getx()");
    return geom().vx_vec();
  }
  vd Duct::getValueT(Varnr v,d timeinst) const{
    TRACE(15,"Duct::getValueT()");
    vd res=getValue(v,0);	// Time-averaged part
    d omg=gc->getomg();
    us Nf=gc->Nf();
    for(us i=1; i<Nf+1 ;i++){
      res+=real(getValueC(v,i)*exp((omg*I)*((d) i)*timeinst));
    }
    return res;
  }
  vd Duct::getValue(Varnr v,us freqnr) const {
    TRACE(15,"Duct::getValue("<<(int)v<<","<<freqnr<<")");
    if(freqnr>=gc->Ns())
      throw MyError("Illegal frequency number");
    const us nCells=geom().nCells();

    vd res(nCells);
    for(us i=0;i<nCells;i++){
      res(i)=cells[i]->getValue(v,freqnr);
    }
    return res;
  }
  vc Duct::getValueC(Varnr v,us freqnr) const {
    TRACE(15,"Duct::getValueC("<<(int)v<<","<<freqnr<<")");
    const us nCells=geom().nCells();
    if(freqnr>gc->Nf() || freqnr<1)
      throw MyError("Illegal frequency number");
    // VARTRACE(15,getNDofs());
    vc res=(1.0+0.0*I)*getValue(v,2*freqnr-1);
    res+=I*getValue(v,2*freqnr);
    return res;
  }
  vd Duct::getErrorAt(us eqnr,us freqnr) const {
    const us& nCells=getNCells();
    vd er(nCells,fillwith::zeros);
    assert(eqnr<getNDofs());
    for(us i=0;i<nCells;i++){
      er(i)=(cells[i]->errorAt(eqnr))(eqnr);
    }
    return er;
  }
  void Duct::resetHarmonics(){
    for(auto v=cells.begin();v!=cells.end();v++){
      (*v)->resetHarmonics();
    }
  }
  void Duct::setPhaseContraint(PhaseConstraint pc){
    TRACE(15,"Duct::setPhaseContraint(PhaseConstraint v)");
    // Sanity checks on Phase constraint:
    if(pc.var == Varnr::p || pc.var == Varnr::m || pc.var==Varnr::T) {
      delete pc_;
      pc_=new PhaseConstraint(pc);
    }
    else {
      WARN("Phase constraint on " << toString(pc.var) << " cannot be set in a Duct");
      throw MyError("Phase constraint cannot be set. See warning messages for details.");
    }
  }
  int Duct::providePhaseDof() const {
    TRACE(20,"Duct::providePhaseDof()");
    if(!pc_) {
      TRACE(20,"No phase contraint in " << getName());
      return -1;
    }
    if(pc_->freqnr>Gc().Ns()-1){
      WARN("Frequency number given for phase constraint is "
	   "too high! Given: " << pc_->freqnr);
      throw MyError("Illegal frequency given. See warnings.");
    }
    if(pc_->var == Varnr::p) 
      return bcCell(pc_->pos).p().getDofNr()+pc_->freqnr;
    if(pc_->var == Varnr::m) 
      return bcCell(pc_->pos).mbc().getDofNr()+pc_->freqnr;
    if(pc_->var == Varnr::T) 
      return bcCell(pc_->pos).Tbc().getDofNr()+pc_->freqnr;
    assert(false);
  }
  d Duct::phaseDofValue() const {
    TRACE(15,"Duct::PhaseDofValue()");
    assert(pc_);
    if(pc_->var == Varnr::p) 
      return bcCell(pc_->pos).p()(pc_->freqnr);
    if(pc_->var == Varnr::m) 
      return bcCell(pc_->pos).mbc()(pc_->freqnr);
    if(pc_->var == Varnr::T) 
      return bcCell(pc_->pos).Tbc()(pc_->freqnr);
    assert(false);
  }
  void Duct::dmtotdx(vd& dmtotdx_) const{
    TRACE(15,"Duct::dmtotdx()");
    us ncell=cells.size(),Neq;
    us rhodof;
    for(auto cell=cells.begin();cell!=cells.end();cell++){
      rhodof=(*cell)->rho().getDofNr();
      dmtotdx_(rhodof)=(*cell)->vVf;
    }
  }
  void Duct::show(us detailnr) const {
    cout << "++++++++++++Duct name: "<< getName() << " ++++++++++++++++\n";
    cout << "Type: " << typeid(*this).name() <<" with number "<<getNumber()<< ".\n";
    cout << "********************************************************************************\n";
    cout << "Geometry: \n";
    assert(cells.size()!=0);
    if(detailnr>2){
      geom().show();
    }
    if(detailnr>3)
      this->showVertices(detailnr);
  }
  void Duct::showVertices(us showvertices) const {
    for(us i=0;i<cells.size();i++)
      cells[i]->show(showvertices);
  }
  void Duct::jac(Jacobian& tofill) const{			// Return Jacobian matrix of error operator
    // sdmat Duct::Jac(){			// Return Jacobian matrix of error operator    
    TRACE(8," Duct::Jac() for Ductment "<< getNumber() << ".");
    const us& Ns=gc->Ns();
    us nCell=cells.size();

    for(us j=0;j<nCell;j++){			   // Fill the Jacobian
      TRACE(3,"Obtaining cell Jacobian...");
      cells.at(j)->jac(tofill);
    }	// end for
	// cout <<"Segment" << getNumber() <<" Jacobian done. Jac is:\n"<< Jacobian;
	// cout << "Number of colums in this jacobian" << Jacobian.n_cols<<"\n";
    TRACE(8,"Ductment Jacobian done.");
  }
  vd Duct::getRes() const {
    TRACE(8,"Duct::GetRes()");
    assert(cells.size()!=0);
    assert(gc!=nullptr);
    vd Result(getNDofs(),fillwith::zeros);
    us nCell=cells.size();    
    us Ns=gc->Ns();
    us vndofs,curpos=0;
    for(us k=0; k<nCell;k++) {
      vndofs=cells.at(k)->getNDofs();
      Result.subvec(curpos,curpos+vndofs-1)=cells[k]->getRes();
      curpos+=vndofs;
    }
    return Result;
  }
  vd Duct::error() const{
    TRACE(8,"Duct::Error()");
    assert(cells.size()!=0);
    assert(gc!=nullptr);
    vd error(getNEqs(),fillwith::zeros);
    us nCell=cells.size();    
    us Ns=gc->Ns();
    us vneqs,curpos=0;

    for(us k=0; k<nCell;k++) {
      vneqs=cells.at(k)->getNEqs();
      error.subvec(curpos,curpos+vneqs-1)=cells[k]->error();
      curpos+=vneqs;
    }
    TRACE(10,"Filling error for seg nr " << getNumber() << " done." );
    return error;
  }
  void Duct::domg(vd& domg) const{
    TRACE(8,"Duct::Error()");
    for(const Cell* cell:cells) {
      cell->domg(domg);
    }
  }

  void Duct::setRes(const vd& res){
    TRACE(8,"Duct::SetRes()");
    assert(res.size()==getNDofs());
    // const us& Neq=(cells[0]).Neq;
    const us& Ns=gc->Ns();
    us celldofs;
    us firstdof=0;
    for(us k=0; k<cells.size();k++) {
      celldofs=cells.at(k)->getNDofs();
      cells.at(k)->setRes(res.subvec(firstdof,firstdof+celldofs-1));
      firstdof+=celldofs;
    }
  }
  // Various set and get methods
  void Duct::setResVar(Varnr v,us i,us freqnr,d value){
    WARN("Func does nothing!");
  }
  void Duct::setResVar(Varnr v,us freqnr,const vd& vals){
    WARN("Func does nothing!");
  }


} /* namespace duct */


