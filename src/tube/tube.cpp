/*
 * lintube.cpp

 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#include "tube.h"
#include "lefttubevertex.h"       // which includes tubebcvertex.h
#include "righttubevertex.h"       // which includes tubebcvertex.h
#include "interpolate.h"
#include "globalconf.h"
#include "geom.h"
#include "exception.h"

// Tried to keep the method definition a bit in order in which a
  // tube is created, including all its components. First a tube is
  // created, which has a geometry and a global
  // configuration. Moreover, the tube has gridpoints, "TubeVertex"
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
    cleanup_vvertex();
  }
  void Tube::cleanup_vvertex(){
    TRACE(25,"Tube::cleanup_vvertex()");
    // for(auto v=vvertex.begin();v!=vvertex.end();v++)
      // delete *v;
    vvertex.clear();
  }

  const TubeVertex& Tube::getTubeVertex(us i) const{
    assert(vvertex.size()>0);
    assert(i<vvertex.size());
    return *static_cast<const TubeVertex*>(vvertex.at(i));
  }
  us Tube::getNCells() const {return geom().nCells();}
  const Geom& Tube::geom() const {return *geom_;}
  void Tube::setDofNrs(us firstdof){
    TRACE(13,"Tube::setDofNrs()");
    assert(vvertex.size()>0);
    for(auto vertex=vvertex.begin();vertex!=vvertex.end();vertex++){
      (*vertex)->setDofNrs(firstdof);
      firstdof+=(*vertex)->getNDofs();
    }
  }
  void Tube::setEqNrs(us firstdof){
    TRACE(13,"Tube::setDofNrs()");
    assert(vvertex.size()>0);
    for(auto vertex=vvertex.begin();vertex!=vvertex.end();vertex++){
      (*vertex)->setEqNrs(firstdof);
      firstdof+=(*vertex)->getNEqs();
    }
  }
  
  us Tube::getNDofs() const {
    TRACE(10,"Tube::getNDofs()");
    us ndofs=0;
    for(auto v=vvertex.begin();v!=vvertex.end();v++)
      ndofs+=(*v)->getNDofs();
    return ndofs;
  }
  us Tube::getNEqs() const {
    TRACE(10,"Tube::getNEqs()");
    us ndofs=0;
    for(auto v=vvertex.begin();v!=vvertex.end();v++)
      ndofs+=(*v)->getNEqs();
    return ndofs;
  }  

  void Tube::init(const tasystem::TaSystem& sys){
    TRACE(13,"Tube::Init()");
    Seg::init(sys);

    cleanup_vvertex();
    TRACE(13,"Filling vertices. Current size:"<<vvertex.size());
      // Left *probable* boundary condition
      // vvertex.emplace_back(new LeftTubeVertex(0,g));
    // WARN("Lot wrong here");
    vvertex.emplace_back(new LeftTubeVertex(0,*this));
    us i;
    for(i=1;i<getNCells()-1;i++)
      vvertex.emplace_back(new TubeVertex(i,*this));
    vvertex.emplace_back(new RightTubeVertex(getNCells()-1,*this));

    us nVertex=vvertex.size();    
    assert(nVertex==getNCells());
    // And initialize again.
    for(i=0;i<vvertex.size();i++){
      TRACE(13,"Starting intialization of Vertex "<< i);
      TubeVertex* thisvertex=vvertex[i];
      TubeVertex* left=NULL;
      TubeVertex* right=NULL;
      if(i<nVertex-1)
        right=vvertex[i+1];
      if(i>0)
        left=vvertex[i-1];
      TRACE(15,"Initializing tube");
      thisvertex->init(left,right);
    } // for
  } // Tube::init(gc)
  d Tube::getRes(us dofnr) const {
    WARN("Not yet implemented!");
  }
  d Tube::getCurrentMass() const{
    TRACE(8,"Tube::getCurrentMass()");
    assert(vvertex.size()>0);
    d mass=0;
    for(auto vertex=vvertex.begin();vertex!=vvertex.end();vertex++){
      TubeVertex& cvertex=*static_cast<TubeVertex*>(*vertex);
      mass+=cvertex.getCurrentMass();
    }
    return mass;
  }
  void Tube::updateNf(){
    TRACE(18,"Tube::updateNf()");
    assert(vvertex.size()>0);
    for(auto v=vvertex.begin();v!=vvertex.end();v++){
      (*v)->updateNf();
    }
  }
  const TubeBcVertex& Tube::leftVertex() const{
    TRACE(3,"Tube::leftVertex()");
    assert(vvertex[0]);
    return static_cast<const TubeBcVertex&>(*vvertex[0]);
  }
  const TubeBcVertex& Tube::rightVertex() const{
    TRACE(3,"Tube::rightVertex()");
    return static_cast<const TubeBcVertex&>(**(vvertex.end()-1));
  }
  vd Tube::getx() const {
    TRACE(10,"Tube::getx()");
    checkInit();
    return geom().vx_vec();
  }
    
  vd Tube::getValue(varnr v,us freqnr) const throw(std::exception) {
    TRACE(10,"Tube::getValue("<<(int)v<<","<<freqnr<<")");
    checkInit();
    if(freqnr>=gc->Ns())
      throw MyError("Illegal frequency number");
    const us nCells=geom().nCells();
    // VARTRACE(15,getNDofs());

    vd res(nCells+2);
    for(us i=0;i<nCells;i++){
      res(i+1)=vvertex[i]->getValue(v,freqnr);
    }
    res(0)=leftVertex().getValueBc(v,freqnr);
    res(nCells+1)=rightVertex().getValueBc(v,freqnr);
    return res;
  }
  vc Tube::getValueC(varnr v,us freqnr) const throw(std::exception) {
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
      er(i)=(vvertex[i]->errorAt(eqnr))(eqnr);
    }
    return er;
  }
  void Tube::resetHarmonics(){
    for(auto v=vvertex.begin();v!=vvertex.end();v++){
      (*v)->resetHarmonics();
    }
  }
  void Tube::dmtotdx(vd& dmtotdx_) const{
    TRACE(15,"Tube::dmtotdx()");
    us nvertex=vvertex.size(),Neq;
    us rhodof;
    for(auto v=vvertex.begin();v!=vvertex.end();v++){
      auto &cvertex=*static_cast<TubeVertex*>(*v);
      rhodof=cvertex.rho().getDofNr();
      dmtotdx_(rhodof)=cvertex.localGeom().vVf;
    }
  }

  vd Tube::Htot() const throw(std::exception){
    TRACE(15,"Tube::Htot()");
    
    us nvertex=vvertex.size();
    vd Htot(nvertex);
    for(us i=0;i<nvertex;i++){
      Htot(i)=vvertex[i]->Htot();
    }
    return Htot;
  }

  void Tube::show(us detailnr) const {
    cout << "++++++++++++Tube name: "<< getName() << " ++++++++++++++++\n";
    cout << "Type: " << getType() <<" with number "<<getNumber()<< ".\n";
    cout << "********************************************************************************\n";
    cout << "Geometry: \n";
    assert(vvertex.size()!=0);
    if(detailnr>2){
      geom().show();
    }
    if(detailnr>3)
      this->showVertices(detailnr);
  }
  void Tube::showVertices(us showvertices) const {
    for(us i=0;i<vvertex.size();i++)
      vvertex[i]->show(showvertices);
  }
  void Tube::jac(Jacobian& tofill) const{			// Return Jacobian matrix of error operator
    // sdmat Tube::Jac(){			// Return Jacobian matrix of error operator    
    TRACE(8," Tube::Jac() for Tubement "<< getNumber() << ".");
    const us& Ns=gc->Ns();
    us nVertex=vvertex.size();

    for(us j=0;j<nVertex;j++){			   // Fill the Jacobian
      TRACE(3,"Obtaining vertex Jacobian...");
      vvertex.at(j)->jac(tofill);
    }	// end for
    // cout <<"Segment" << getNumber() <<" Jacobian done. Jac is:\n"<< Jacobian;
    // cout << "Number of colums in this jacobian" << Jacobian.n_cols<<"\n";
    TRACE(8,"Tubement Jacobian done.");
  }
  vd Tube::getRes() const {
    TRACE(8,"Tube::GetRes()");
    assert(vvertex.size()!=0);
    assert(gc!=NULL);
    vd Result(getNDofs(),fillwith::zeros);
    us nVertex=vvertex.size();    
    us Ns=gc->Ns();
    us vndofs,curpos=0;
    for(us k=0; k<nVertex;k++) {
      vndofs=vvertex.at(k)->getNDofs();
      Result.subvec(curpos,curpos+vndofs-1)=vvertex[k]->getRes();
      curpos+=vndofs;
    }
    return Result;
  }
  vd Tube::error() const{
    TRACE(8,"Tube::Error()");
    assert(vvertex.size()!=0);
    assert(gc!=NULL);
    vd error(getNEqs(),fillwith::zeros);
    us nVertex=vvertex.size();    
    us Ns=gc->Ns();
    us vneqs,curpos=0;

    for(us k=0; k<nVertex;k++) {
      vneqs=vvertex.at(k)->getNEqs();
      error.subvec(curpos,curpos+vneqs-1)=vvertex[k]->error();
      curpos+=vneqs;
    }
    TRACE(10,"Filling error for seg nr " << getNumber() << " done." );
    return error;
  }
  void Tube::domg(vd& domg_v) const{
    TRACE(8,"Tube::Error()");
    const us& Ns=gc->Ns();
    us nVertex=vvertex.size();    
    // vd domg(getNDofs(),fillwith::zeros);
    us vndofs,curpos=0;
    for(us k=0; k<nVertex;k++) {
      vvertex[k]->domg(domg_v);
    }
  }

  void Tube::setRes(const vd& res){
    TRACE(8,"Tube::SetRes()");
    assert(res.size()==getNDofs());
    // const us& Neq=(vvertex[0]).Neq;
    const us& Ns=gc->Ns();
    us vertexdofs;
    us firstdof=0;
    for(us k=0; k<vvertex.size();k++) {
      vertexdofs=vvertex.at(k)->getNDofs();
      vvertex.at(k)->setRes(res.subvec(firstdof,firstdof+vertexdofs-1));
      firstdof+=vertexdofs;
    }
  }

  // Various set and get methods
  void Tube::setResVar(varnr v,us i,us freqnr,d value){
    WARN("Func does nothing!");
  }
  void Tube::setResVar(varnr v,us freqnr,const vd& vals){
    TRACE(15,"Tube::setResVar()");
    if(v==varnr::p){
      assert(vals.size()==geom().nCells()+1);
    }
    else{
      assert(vals.size()==geom().nCells()+2);
    }
      
    for(auto v=vvertex.begin();v!=vvertex.end();v++){

    }
  }

  void Tube::setRes(const Seg& otherseg){
    TRACE(20,"Tube::setRes(othertube)");
    const Tube& other=asTube_const(otherseg);
    // Sanity checks
    assert(vvertex.size()!=0);
    // Necessary to let it work
    assert(gc->Ns()==other.gc->Ns());

  
    for(auto v=vvertex.begin();v!=vvertex.end();v++){
      TubeVertex& thisvertex=*static_cast<TubeVertex*>(*v);
      d vx=thisvertex.localGeom().vx;
      d xL=thisvertex.localGeom().xL;
      thisvertex.setResVar(varnr::rho,other.interpolateResMid(varnr::rho,vx));
      thisvertex.setResVar(varnr::U,other.interpolateResMid(varnr::U,vx));
      thisvertex.setResVar(varnr::T,other.interpolateResMid(varnr::T,vx));
      thisvertex.setResVar(varnr::Ts,other.interpolateResMid(varnr::Ts,vx));
      thisvertex.setResVar(varnr::p,other.interpolateResStaggered(varnr::p,xL));
      WARN("boundaries todo");
      // if(v==(vvertex.end()-1)){
      //   if(bcRight){
      //     TubeVertex& othervertex=*static_cast<TubeVertex*>(*(other.vvertex.end()-1));          
      //     thisvertex.setpR(othervertex.pR());
      //     TRACE(5,"Copying pR");          
      //   }
      // }
    } // for

  }
  vd Tube::interpolateResStaggered(varnr v,d x) const{
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
    // vd left=leftvertex.getRes(v)();
    // vd right=rightvertex.getRes(v)();
    // d xleft=leftvertex.localGeom().xL;
    // d xright=rightvertex.localGeom().xL;
    // d relpos=(x-xleft)/(xright-xleft);
    // VARTRACE(5,relpos);
    // return math_common::linearInterpolate(left,right,relpos);
  }

  vd Tube::interpolateResMid(varnr v,d x) const{
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
    const TubeVertex& leftvertex=*vvertex[ileft];
    const TubeVertex& rightvertex=*vvertex[iright];
    // vd left=leftvertex.getRes(v)();
    // vd right=rightvertex.getRes(v)();
    // d xleft=leftvertex.localGeom().vx;
    // d xright=rightvertex.localGeom().vx;
    // d relpos=(x-xleft)/(xright-xleft);
    // VARTRACE(5,relpos);
    // return math_common::linearInterpolate(left,right,relpos);
  }

} /* namespace tube */


