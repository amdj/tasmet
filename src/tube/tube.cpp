/*
 * lintube.cpp

 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#include "tube.h"
#include "leftcell.h"       // which includes tubebccell.h
#include "rightcell.h"       // which includes tubebccell.h
#include "globalconf.h"
#include "geom.h"
#include "exception.h"
#include "utils.h"

// Tried to keep the method definition a bit in order in which a
  // tube is created, including all its components. First a tube is
  // created, which has a geometry and a global
  // configuration. Moreover, the tube has gridpoints, "Cell"
  // instants. Of these, a tube has gp of them, stored in a vector. In
  // each gridpoint, variables live, which represent the current
  // solution. Moreover, we have equations in each gridpoint. More
  // precisely, in the final solution the continuity, momentum, energy
  // and a suitable equation of state should hold.

namespace tube {
  using tasystem::PhaseConstraint;
  using tasystem::TaSystem;
  using tasystem::Globalconf;
  using tasystem::Jacobian;


  Tube::Tube(const Geom& geom)
    :Seg(),geom_(geom.copy())
  {
    TRACE(13,"Tube constructor()...");
  }
  Tube::Tube(const Tube& other,const TaSystem& sys):
    Seg(other,sys),
    geom_(other.geom().copy()){
    if(other.pc_)
      this->pc_=new PhaseConstraint(*other.pc_);
    TRACE(13,"Filling vertices. Current size:"<<cells.size());
      // Left *probable* boundary condition
      // cells.emplace_back(new LeftCell(0,g));
    // WARN("Lot wrong here");
    cells.emplace_back(new LeftCell(0,*this));
    us i;
    for(i=1;i<getNCells()-1;i++)
      cells.emplace_back(new Cell(i,*this));
    cells.emplace_back(new RightCell(getNCells()-1,*this));

    us nCell=cells.size();    
    assert(nCell==getNCells());
    // And initialize again.
    for(i=0;i<cells.size();i++){
      TRACE(13,"Starting intialization of Cell "<< i);
      Cell* thiscell=cells[i];
      Cell* left=nullptr;
      Cell* right=nullptr;
      if(i<nCell-1)
        right=cells[i+1];
      if(i>0)
        left=cells[i-1];
      TRACE(15,"Initializing tube");
      thiscell->init(left,right);
    } // for
  } // Tube::Tube(copy)

  Tube::~Tube(){
    TRACE(25,"~Tube()");
    delete geom_;
    delete pc_;
    utils::purge(cells);
  }
  const Cell& Tube::getCell(int i) const{
    if(i>=0)
      return *cells.at(i);
    else
      return *cells.at(cells.size()+i);
  }
  const Cell& Tube::operator[](us i) const{
    return *cells[i];
  }
  us Tube::getNCells() const {return geom().nCells();}
  const Geom& Tube::geom() const {return *geom_;}
  void Tube::setDofNrs(us firstdof){
    TRACE(13,"Tube::setDofNrs("<<firstdof<<")");
    assert(cells.size()>0);
    for(Cell* cell: cells){
      cell->setDofNrs(firstdof);
      firstdof+=cell->getNDofs();
    }
  }
  void Tube::setEqNrs(us firsteq){
    TRACE(13,"Tube::SetEqNrs("<<firsteq<<")");
    assert(cells.size()>0);
    for(auto cell=cells.begin();cell!=cells.end();cell++){
      (*cell)->setEqNrs(firsteq);
      firsteq+=(*cell)->getNEqs();
    }
  }
  
  us Tube::getNDofs() const {
    TRACE(10,"Tube::getNDofs()");
    us ndofs=0;
    for(auto v=cells.begin();v!=cells.end();v++)
      ndofs+=(*v)->getNDofs();
    return ndofs;
  }
  us Tube::getNEqs() const {
    TRACE(10,"Tube::getNEqs()");
    us ndofs=0;
    for(auto v=cells.begin();v!=cells.end();v++)
      ndofs+=(*v)->getNEqs();
    return ndofs;
  }  

  d Tube::getMass() const{
    TRACE(8,"Tube::getMass()");
    assert(cells.size()>0);
    d mass=0;
    for(auto cell=cells.begin();cell!=cells.end();cell++){
      mass+=(*cell)->getMass();
    }
    return mass;
  }
  void Tube::updateNf(){
    TRACE(18,"Tube::updateNf()");
    assert(cells.size()>0);
    for(auto v: cells){
      v->updateNf();
    }
  }
  // vc Tube::heatQ() const {
  //   TRACE(15,"Tube::heatQ()");
  //   vc res(getNCells());
  //   for(us i=0; i< getNCells(); i++) {
  //     res(i)=static_cast<const HopkinsHeatSource&>(getHeatSource()).HeatTransferCoefQ(*cells[i])(1);
  //   }
  //   return res;
  // }
  // vc Tube::heatH() const {
  //   TRACE(15,"Tube::heatQ()");
  //   vc res(getNCells());
  //   for(us i=0; i< getNCells(); i++) {
  //     res(i)=static_cast<const HopkinsHeatSource&>(getHeatSource()).HeatTransferCoefH(*cells[i])(1);
  //   }
  //   return res;
  // }
  const BcCell& Tube::bcCell(Pos p) const{
    TRACE(3,"Tube::bcCell()");
    if(p==Pos::left){
      assert(cells[0]);
      return static_cast<const BcCell&>(*cells[0]);
    }
    else{
      return static_cast<const BcCell&>(**(cells.end()-1));
    }
  }
  vd Tube::getx() const {
    TRACE(10,"Tube::getx()");
    checkInit();
    return geom().vx_vec();
  }
  vd Tube::getValueT(Varnr v,d timeinst) const{
    TRACE(15,"Tube::getValueT()");
    vd res=getValue(v,0);	// Time-averaged part
    // omega*T=2*pi
    us Nf=gc->Nf();
    for(us i=1; i<Nf+1 ;i++){
      res+=real(getValueC(v,i)*exp((2.0*I*number_pi)*((d) i)*timeinst));
    }
    return res;
  }
  vd Tube::getValue(Varnr v,us freqnr) const {
    TRACE(15,"Tube::getValue("<<(int)v<<","<<freqnr<<")");
    checkInit();
    if(freqnr>=gc->Ns())
      throw MyError("Illegal frequency number");
    const us nCells=geom().nCells();

    vd res(nCells);
    for(us i=0;i<nCells;i++){
      res(i)=cells[i]->getValue(v,freqnr);
    }
    return res;
  }
  vc Tube::getValueC(Varnr v,us freqnr) const {
    TRACE(15,"Tube::getValueC("<<(int)v<<","<<freqnr<<")");
    const us nCells=geom().nCells();
    if(freqnr>gc->Nf() || freqnr<1)
      throw MyError("Illegal frequency number");
    // VARTRACE(15,getNDofs());
    vc res=(1.0+0.0*I)*getValue(v,2*freqnr-1);
    res+=I*getValue(v,2*freqnr);
    return res;
  }
  vd Tube::getErrorAt(us eqnr,us freqnr) const {
    const us& nCells=getNCells();
    vd er(nCells,fillwith::zeros);
    assert(eqnr<getNDofs());
    for(us i=0;i<nCells;i++){
      er(i)=(cells[i]->errorAt(eqnr))(eqnr);
    }
    return er;
  }
  void Tube::resetHarmonics(){
    for(auto v=cells.begin();v!=cells.end();v++){
      (*v)->resetHarmonics();
    }
  }
  void Tube::setPhaseContraint(PhaseConstraint pc){
    TRACE(15,"Tube::setPhaseContraint(PhaseConstraint v)");
    // Sanity checks on Phase constraint:
    if(pc.var == Varnr::p || pc.var == Varnr::m || pc.var==Varnr::T) {
      delete pc_;
      pc_=new PhaseConstraint(pc);
    }
    else {
      WARN("Phase constraint on " << toString(pc.var) << " cannot be set in a Tube");
      throw MyError("Phase constraint cannot be set. See warning messages for details.");
    }
  }
  int Tube::providePhaseDof() const {
    TRACE(20,"Tube::providePhaseDof()");
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
  d Tube::phaseDofValue() const {
    TRACE(15,"Tube::PhaseDofValue()");
    assert(pc_);
    if(pc_->var == Varnr::p) 
      return bcCell(pc_->pos).p()(pc_->freqnr);
    if(pc_->var == Varnr::m) 
      return bcCell(pc_->pos).mbc()(pc_->freqnr);
    if(pc_->var == Varnr::T) 
      return bcCell(pc_->pos).Tbc()(pc_->freqnr);
    assert(false);
  }
  void Tube::dmtotdx(vd& dmtotdx_) const{
    TRACE(15,"Tube::dmtotdx()");
    us ncell=cells.size(),Neq;
    us rhodof;
    for(auto cell=cells.begin();cell!=cells.end();cell++){
      rhodof=(*cell)->rho().getDofNr();
      dmtotdx_(rhodof)=(*cell)->vVf;
    }
  }
  void Tube::show(us detailnr) const {
    cout << "++++++++++++Tube name: "<< getName() << " ++++++++++++++++\n";
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
  void Tube::showVertices(us showvertices) const {
    for(us i=0;i<cells.size();i++)
      cells[i]->show(showvertices);
  }
  void Tube::jac(Jacobian& tofill) const{			// Return Jacobian matrix of error operator
    // sdmat Tube::Jac(){			// Return Jacobian matrix of error operator    
    TRACE(8," Tube::Jac() for Tubement "<< getNumber() << ".");
    const us& Ns=gc->Ns();
    us nCell=cells.size();

    for(us j=0;j<nCell;j++){			   // Fill the Jacobian
      TRACE(3,"Obtaining cell Jacobian...");
      cells.at(j)->jac(tofill);
    }	// end for
    // cout <<"Segment" << getNumber() <<" Jacobian done. Jac is:\n"<< Jacobian;
    // cout << "Number of colums in this jacobian" << Jacobian.n_cols<<"\n";
    TRACE(8,"Tubement Jacobian done.");
  }
  vd Tube::getRes() const {
    TRACE(8,"Tube::GetRes()");
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
  vd Tube::error() const{
    TRACE(8,"Tube::Error()");
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
  void Tube::domg(vd& domg) const{
    TRACE(8,"Tube::Error()");
    for(const Cell* cell:cells) {
      cell->domg(domg);
    }
  }

  void Tube::setRes(const vd& res){
    TRACE(8,"Tube::SetRes()");
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
  void Tube::setResVar(Varnr v,us i,us freqnr,d value){
    WARN("Func does nothing!");
  }
  void Tube::setResVar(Varnr v,us freqnr,const vd& vals){
    WARN("Func does nothing!");
  }


} /* namespace tube */


