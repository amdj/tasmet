/*
 * lintube.cpp

 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#include "tube.h"
#include "tubebcvertex.h"
#include "interpolate.h"
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

  Tube::Tube(const Geom& geom):Seg(geom){
    // Fill vector of gridpoints with data:
    TRACE(13,"Tube constructor()...");
  }
  void Tube::cleanup(){
    TRACE(13,"Tube::cleanup()");
    bcLeft.reset();
    bcRight.reset();
    for(auto v=vvertex.begin();v!=vvertex.end();v++)
      delete *v;
    vvertex.clear();
  }
  Tube::Tube(const Tube& other):Seg(other){
    copyTube(other);
  }
  Tube& Tube::operator=(const Tube& other){
    TRACE(13,"Tube copy assignment");
    cleanup();
    Seg::operator=(other);
    copyTube(other);
    return *this;
  }
  const TubeVertex& Tube::getTubeVertex(us i) const{
    assert(vvertex.size()>0);
    assert(i<vvertex.size());
    return *static_cast<const TubeVertex*>(vvertex.at(i));
  }
  void Tube::show(us showvertices) const {
    TRACE(18,"Tube::show()");
    cout << "++++++++++++Tube name: "<< getName() << " ++++++++++++++++\n";
    if(bcLeft){
      cout << "Left side contains internal boundary condition of type " << bcLeft->getType() << ".\n";
    }
    if(bcRight){
      cout << "Right side contains internal boundary condition of type " << bcRight->getType() << ".\n";
    }
    Seg::show(showvertices);
    cout << "********************************************************************************\n";
  }
  vd Tube::dragCoefVec(us i) const {
    return zeros<vd>(geom().nCells());
  }
  void Tube::copyTube(const Tube& other){
    TRACE(13,"Tube::copyTube()");
    cleanup();
    if(other.bcLeft){
      TRACE(12,"Other has bc connected left side");
      bcLeft.reset(other.bcLeft->copy());
    }
    if(other.bcRight)
      bcRight.reset(other.bcRight->copy());
  }
  void Tube::setDofNrs(us firstdof){
    TRACE(13,"Tube::setDofNrs()");
    assert(vvertex.size()>0);
    for(auto vertex=vvertex.begin();vertex!=vvertex.end();vertex++){
      TubeVertex& v=*static_cast<TubeVertex*>(*vertex);
      v.setDofNrs(firstdof);
      firstdof+=v.getNDofs();
    }
  }
  void Tube::setEqNrs(us firstdof){
    TRACE(13,"Tube::setDofNrs()");
    assert(vvertex.size()>0);
    for(auto vertex=vvertex.begin();vertex!=vvertex.end();vertex++){
      TubeVertex& v=*static_cast<TubeVertex*>(*vertex);
      v.setEqNrs(firstdof);
      firstdof+=v.getNEqs();
    }
  }
  
  void Tube::addBc(const TubeBcVertex& bc){
    TRACE(14,"Tube::addBc(bc)");
    if(bc.connectPos()==connectpos::left)
      {
	TRACE(12,"Bc connected left side");
	bcLeft.reset(bc.copy());
	if(bcLeft)
	  TRACE(12,"bcLeft is now "<< bool(bcLeft));
      }
    else if(bc.connectPos()==connectpos::right)
      {
	TRACE(12,"Bc connected right side");
	bcRight.reset(bc.copy());
      }
    else
      {      
	cout << "WARNING: bconnectbc(): Bc  not understood!\n";
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
  TubeVertex* Tube::leftTubeVertex() const{
    TRACE(13,"Tube::leftTubeVertex()");
    if(bcLeft){
      TRACE(13,"Tube::leftTubeVertex() returning a boundary vertex");
      return static_cast<TubeVertex*>(bcLeft->copy());
    }
    else{
      TRACE(13,"Tube::leftTubeVertex() returning an ordinary vertex");
      return new TubeVertex();
    }
  }
  TubeVertex* Tube::rightTubeVertex() const{
    if(bcRight)
      return static_cast<TubeVertex*>(bcRight->copy());
    else
      return new TubeVertex();
  }


  void Tube::init(const tasystem::Globalconf& g){
    TRACE(13,"Tube::Init()");
    Seg::init(g);
    if(vvertex.size()==0){
      TRACE(13,"Filling vertices. Current size:"<<vvertex.size());
      // Left *probable* boundary condition
      vvertex.emplace_back(leftTubeVertex());
      for(us i=1;i<geom().nCells()-1;i++)
    	vvertex.emplace_back(new TubeVertex());
      // Right *probable* boundary condition
      vvertex.emplace_back(rightTubeVertex());

      us nVertex=vvertex.size();    
      assert(nVertex==geom().nCells());

      for(us i=0;i<vvertex.size();i++){
        vvertex.at(i)->initVertex(i,*this);
      }
      // And initialize again.
      for(us i=0;i<vvertex.size();i++){
        TubeVertex* cvertex=static_cast<TubeVertex*>(vvertex[i]);
        TRACE(13,"Starting intialization of Vertex "<< i);
        if(i<nVertex-1) cvertex->setRight(*vvertex[i+1]);
        if(i>0) cvertex->setLeft(*vvertex[i-1]);
        cvertex->initTubeVertex(i,*this);
      } // for
    } // vertex.size==0
    else{
      TRACE(13,"Tube already initialized!");
    }
  } // Tube::init(gc)
  d Tube::getCurrentMass() const{
    TRACE(8,"Tube::getCurrentMass()");
    assert(vvertex.size()>0);
    d mass=0;
    for(auto vertex=vvertex.begin();vertex!=vvertex.end();vertex++){
      TubeVertex& cvertex=*static_cast<TubeVertex*>(*vertex);
      mass+=cvertex.rho(0)*cvertex.lg.vVf;
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
  void Tube::setRes(const SegBase& otherseg){
    TRACE(20,"Tube::setRes(othertube)");
    const Tube& other=asTube_const(otherseg);
    // Sanity checks
    assert(vvertex.size()!=0);
    // Necessary to let it work
    assert(gc->Ns()==other.gc->Ns());

  
    for(auto v=vvertex.begin();v!=vvertex.end();v++){
      TubeVertex& thisvertex=*static_cast<TubeVertex*>(*v);
      d vx=thisvertex.lg.vx;
      d xL=thisvertex.lg.xL;
      thisvertex.rho.set(other.interpolateResMid(varnr::rho,vx));
      thisvertex.U.set(other.interpolateResMid(varnr::U,vx));
      thisvertex.T.set(other.interpolateResMid(varnr::T,vx));
      thisvertex.Ts.set(other.interpolateResMid(varnr::Ts,vx));
      thisvertex.p.set(other.interpolateResStaggered(varnr::p,xL));

      if(v==(vvertex.end()-1)){
        if(bcRight){
          TubeVertex& othervertex=*static_cast<TubeVertex*>(*(other.vvertex.end()-1));          
          thisvertex.setpR(othervertex.pR());
          TRACE(5,"Copying pR");          
        }
      }
    } // for

  }
  vd Tube::interpolateResStaggered(varnr v,d x) const{
    TRACE(2,"Tube::interpolateResStaggered("<<v<<","<<x<<")");
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
    const TubeVertex& leftvertex=*static_cast<TubeVertex*>(vvertex[ileft]);
    const TubeVertex& rightvertex=*static_cast<TubeVertex*>(vvertex[iright]);
    vd left=leftvertex.getRes(v)();
    vd right=rightvertex.getRes(v)();
    d xleft=leftvertex.lg.xL;
    d xright=rightvertex.lg.xL;
    d relpos=(x-xleft)/(xright-xleft);
    VARTRACE(5,relpos);
    return math_common::linearInterpolate(left,right,relpos);
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
    const TubeVertex& leftvertex=*static_cast<TubeVertex*>(vvertex[ileft]);
    const TubeVertex& rightvertex=*static_cast<TubeVertex*>(vvertex[iright]);
    vd left=leftvertex.getRes(v)();
    vd right=rightvertex.getRes(v)();
    d xleft=leftvertex.lg.vx;
    d xright=rightvertex.lg.vx;
    d relpos=(x-xleft)/(xright-xleft);
    VARTRACE(5,relpos);
    return math_common::linearInterpolate(left,right,relpos);
  }
  vd Tube::getResAt(us varNr,us freqnr) const{
    assert(varNr<NVARS);
    return getResAt((varnr) varNr,freqnr);
  }
  vd Tube::getResAt(varnr v,us freqnr) const{
    TRACE(10,"Tube::getResAt("<<(int)v<<","<<freqnr<<")");
    const us nCells=geom().nCells();
    // VARTRACE(15,getNDofs());
    vd res(nCells);
    for(us i=0;i<nCells;i++){
      TubeVertex& cvertex=*static_cast<TubeVertex*>(vvertex[i]);
      res(i)=cvertex.getRes(v,freqnr);
    }
    return res;
  }
  vd Tube::getErrorAt(us eqnr,us freqnr) const{
    const us& nCells=getNCells();
    vd er(nCells,fillwith::zeros);
    assert(eqnr<getNDofs());
    for(us i=0;i<nCells;i++){
      TubeVertex& cvertex=*static_cast<TubeVertex*>(vvertex[i]);
      er(i)=(cvertex.eqs.at(eqnr)->error(cvertex))(freqnr);
    }
    return er;
  }
  void Tube::resetHarmonics(){
    for(auto v=vvertex.begin();v!=vvertex.end();v++){
      auto &cvertex=*static_cast<TubeVertex*>(*v);
      cvertex.resetHarmonics();
    }
  }
  void Tube::dmtotdx(vd& dmtotdx_) const{
    TRACE(15,"Tube::dmtotdx()");
    us nvertex=vvertex.size(),Neq;
    us rhodof;
    for(auto v=vvertex.begin();v!=vvertex.end();v++){
      auto &cvertex=*static_cast<TubeVertex*>(*v);
      rhodof=cvertex.rho.getDofNr();
      dmtotdx_(rhodof)=cvertex.lg.vVf;
    }
  }

  vd Tube::Htot() const{
    TRACE(15,"Tube::Htot()");
    
    us nvertex=vvertex.size();
    vd Htot(nvertex);
    for(us i=0;i<nvertex;i++){
      auto &cvertex=*static_cast<TubeVertex*>(vvertex.at(i));
      Htot(i)=cvertex.e.Htot(cvertex);
    }
    return Htot;
  }
  
  Tube::~Tube(){
    TRACE(15,"~Tube()");
    cleanup();
  }
  

  
} /* namespace tube */


