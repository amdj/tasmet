/*
 * lintube.cpp

 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#include "tube.h"
#include "leftcell.h"       // which includes tubebccell.h
#include "rightcell.h"       // which includes tubebccell.h
#include "interpolate.h"
#include "globalconf.h"
#include "geom.h"
#include "exception.h"

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
  using tasystem::TaSystem;
  using tasystem::Globalconf;
  using tasystem::Jacobian;


  Tube::Tube(const Geom& geom)
    :Seg(),geom_(geom.copy())
  {
    TRACE(13,"Tube constructor()...");
  }
  Tube::Tube(const Tube& other):
    Seg(other),
    geom_(other.geom().copy()){

  }
  Tube::~Tube(){
    TRACE(25,"~Tube()");
    delete geom_;
    cleanup_cells();
  }
  void Tube::cleanup_cells(){
    TRACE(25,"Tube::cleanup_cells()");
    for(auto v=cells.begin();v!=cells.end();v++)
      delete *v;
    cells.clear();
  }
  const Cell& Tube::getCell(us i) const{
    assert(cells.size()>0);
    assert(i<cells.size());
    return *static_cast<const Cell*>(cells.at(i));
  }
  us Tube::getNCells() const {return geom().nCells();}
  const Geom& Tube::geom() const {return *geom_;}
  void Tube::setDofNrs(us firstdof){
    TRACE(13,"Tube::setDofNrs()");
    assert(cells.size()>0);
    for(auto cell=cells.begin();cell!=cells.end();cell++){
      (*cell)->setDofNrs(firstdof);
      firstdof+=(*cell)->getNDofs();
    }
  }
  void Tube::setEqNrs(us firstdof){
    TRACE(13,"Tube::setDofNrs()");
    assert(cells.size()>0);
    for(auto cell=cells.begin();cell!=cells.end();cell++){
      (*cell)->setEqNrs(firstdof);
      firstdof+=(*cell)->getNEqs();
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

  void Tube::init(const tasystem::TaSystem& sys){
    TRACE(13,"Tube::Init()");
    Seg::init(sys);

    cleanup_cells();
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
  } // Tube::init(gc)
  d Tube::getRes(us dofnr) const {
    WARN("Not yet implemented!");
  }
  d Tube::getCurrentMass() const{
    TRACE(8,"Tube::getCurrentMass()");
    assert(cells.size()>0);
    d mass=0;
    for(auto cell=cells.begin();cell!=cells.end();cell++){
      mass+=(*cell)->getCurrentMass();
    }
    return mass;
  }
  void Tube::updateNf(){
    TRACE(18,"Tube::updateNf()");
    assert(cells.size()>0);
    for(auto v=cells.begin();v!=cells.end();v++){
      (*v)->updateNf();
    }
  }
  const BcCell& Tube::leftCell() const{
    TRACE(3,"Tube::leftCell()");
    assert(cells[0]);
    return static_cast<const BcCell&>(*cells[0]);
  }
  const BcCell& Tube::rightCell() const{
    TRACE(3,"Tube::rightCell()");
    return static_cast<const BcCell&>(**(cells.end()-1));
  }
  vd Tube::getx() const {
    TRACE(10,"Tube::getx()");
    checkInit();
    return geom().vx_vec();
  }
    
  vd Tube::getValue(Varnr v,us freqnr) const throw(std::exception) {
    TRACE(10,"Tube::getValue("<<(int)v<<","<<freqnr<<")");
    checkInit();
    if(freqnr>=gc->Ns())
      throw MyError("Illegal frequency number");
    const us nCells=geom().nCells();
    // VARTRACE(15,getNDofs());

    vd res(nCells+2);
    for(us i=0;i<nCells;i++){
      res(i+1)=cells[i]->getValue(v,freqnr);
    }
    res(0)=leftCell().getValueBc(v,freqnr);
    res(nCells+1)=rightCell().getValueBc(v,freqnr);
    return res;
  }
  vc Tube::getValueC(Varnr v,us freqnr) const throw(std::exception) {
    TRACE(10,"Tube::getResAt("<<(int)v<<","<<freqnr<<")");
    const us nCells=geom().nCells();
    if(freqnr>gc->Nf() || freqnr<1)
      throw MyError("Illegal frequency number");
    // VARTRACE(15,getNDofs());
    vc res(nCells+2);
    res=getValue(v,2*freqnr-1)+I*getValue(v,2*freqnr);
    return res;
  }
  vd Tube::getErrorAt(us eqnr,us freqnr) const throw(std::exception){
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
  void Tube::dmtotdx(vd& dmtotdx_) const{
    TRACE(15,"Tube::dmtotdx()");
    us ncell=cells.size(),Neq;
    us rhodof;
    for(auto cell=cells.begin();cell!=cells.end();cell++){
      rhodof=(*cell)->rho().getDofNr();
      dmtotdx_(rhodof)=(*cell)->vVf;
    }
  }

  // vd Tube::Htot() const throw(std::exception){
  //   TRACE(15,"Tube::Htot()");
    
  //   us ncell=cells.size();
  //   vd Htot(ncell);
  //   for(us i=0;i<ncell;i++){
  //     Htot(i)=cells[i]->Htot();
  //   }
  //   return Htot;
  // }

  void Tube::show(us detailnr) const {
    cout << "++++++++++++Tube name: "<< getName() << " ++++++++++++++++\n";
    cout << "Type: " << getType() <<" with number "<<getNumber()<< ".\n";
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
  void Tube::domg(vd& domg_v) const{
    TRACE(8,"Tube::Error()");
    const us& Ns=gc->Ns();
    us nCell=cells.size();    
    // vd domg(getNDofs(),fillwith::zeros);
    us vndofs,curpos=0;
    for(us k=0; k<nCell;k++) {
      cells[k]->domg(domg_v);
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
    TRACE(15,"Tube::setResVar()");
    if(v==Varnr::p){
      assert(vals.size()==geom().nCells()+1);
    }
    else{
      assert(vals.size()==geom().nCells()+2);
    }
      
    for(auto v=cells.begin();v!=cells.end();v++){

    }
  }

  void Tube::setRes(const Seg& otherseg){
    TRACE(20,"Tube::setRes(othertube)");
    const Tube& other=asTube_const(otherseg);
    // Sanity checks
    assert(cells.size()!=0);
    // Necessary to let it work
    assert(gc->Ns()==other.gc->Ns());

  
    for(auto v=cells.begin();v!=cells.end();v++){
      Cell& thiscell=*(*v);
      d vx=thiscell.vx;
      d xL=thiscell.xL;
      thiscell.setResVar(Varnr::rho,other.interpolateResMid(Varnr::rho,vx));
      thiscell.setResVar(Varnr::U,other.interpolateResMid(Varnr::U,vx));
      thiscell.setResVar(Varnr::T,other.interpolateResMid(Varnr::T,vx));
      thiscell.setResVar(Varnr::Ts,other.interpolateResMid(Varnr::Ts,vx));
      thiscell.setResVar(Varnr::p,other.interpolateResStaggered(Varnr::p,xL));
      WARN("boundaries todo");
      // if(v==(cells.end()-1)){
      //   if(bcRight){
      //     Cell& othercell=*static_cast<Cell*>(*(other.cells.end()-1));          
      //     thiscell.setpR(othercell.pR());
      //     TRACE(5,"Copying pR");          
      //   }
      // }
    } // for

  }
  vd Tube::interpolateResStaggered(Varnr v,d x) const{
    TRACE(2,"Tube::interpolateResStaggered("<<v<<","<<x<<")");
    WARN("out of order!");
    us leftpos=0;
    assert(x>=0);
    us iright=0,ileft;
    while(geom().x(iright)<=x && iright<geom().nCells()-1)
      iright++;
    VARTRACE(5,iright);
    if(iright>0)
      ileft=iright-1;
    else{
      ileft=0;
      iright=1;
    }
    VARTRACE(2,ileft);
    VARTRACE(2,iright);
    // vd left=leftcell.getRes(v)();
    // vd right=rightcell.getRes(v)();
    // d xleft=leftcell.xL;
    // d xright=rightcell.xL;
    // d relpos=(x-xleft)/(xright-xleft);
    // VARTRACE(5,relpos);
    // return math_common::linearInterpolate(left,right,relpos);
  }

  vd Tube::interpolateResMid(Varnr v,d x) const{
    TRACE(2,"Tube::interpolateResMid("<<v<<","<<x<<")");
    us leftpos=0;
    assert(x>=0);
    us iright=0,ileft;
    while(geom().vx(iright)<=x && iright<geom().nCells()-1)
      iright++;
    VARTRACE(5,iright);
    if(iright>0)
      ileft=iright-1;
    else{
      ileft=0;
      iright=1;
    }
    VARTRACE(2,ileft);
    VARTRACE(2,iright);
    const Cell& leftcell=*cells[ileft];
    const Cell& rightcell=*cells[iright];
    // vd left=leftcell.getRes(v)();
    // vd right=rightcell.getRes(v)();
    // d xleft=leftcell.vx;
    // d xright=rightcell.vx;
    // d relpos=(x-xleft)/(xright-xleft);
    // VARTRACE(5,relpos);
    // return math_common::linearInterpolate(left,right,relpos);
  }

} /* namespace tube */


