/*
 * lintube.cpp

 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#include "tube.h"
#include "tubebcvertex.h"

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
    vvertex.clear();
  }
  void Tube::setRes(const SegBase& otherseg){
    TRACE(20,"Tube::setRes(othertube)");
    const Tube& other=asTube_const(otherseg);
    // Sanity checks
    assert(vvertex.size()!=0);
    // Necessary to let it work
    assert(vvertex.size()==other.vvertex.size());
    assert(gc->Ns==other.gc->Ns);
    auto otherv=other.vvertex.begin();

    for(auto v=vvertex.begin();v!=vvertex.end();v++){
      TubeVertex& thisvertex=*static_cast<TubeVertex*>(v->get());
      TubeVertex& othervertex=*static_cast<TubeVertex*>(otherv->get());
      thisvertex.rho=othervertex.rho;
      // TRACE(15,"Other rho:"<< othervertex.rho(0));
      thisvertex.U=othervertex.U;
      thisvertex.T=othervertex.T;      
      thisvertex.Ts=othervertex.Ts;
      if(thisvertex.p.getDofNr()!=-1) // Then its invalid
        thisvertex.p=othervertex.pL();
      if(v==(vvertex.end()-1)){
        if(bcRight){
          thisvertex.setpR(othervertex.pR());
          TRACE(25,"Copying pR");          
        }
      }
      otherv++;
    } // for
    
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
      TubeVertex& v=*static_cast<TubeVertex*>(vertex->get());
      v.setDofNrs(firstdof);
      firstdof+=v.getNDofs();
    }
  }
  void Tube::setEqNrs(us firstdof){
    TRACE(13,"Tube::setDofNrs()");
    assert(vvertex.size()>0);
    for(auto vertex=vvertex.begin();vertex!=vvertex.end();vertex++){
      TubeVertex& v=*static_cast<TubeVertex*>(vertex->get());
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
      ndofs+=(v->get())->getNDofs();
    return ndofs;
  }
  us Tube::getNEqs() const {
    TRACE(10,"Tube::getNEqs()");
    us ndofs=0;
    for(auto v=vvertex.begin();v!=vvertex.end();v++)
      ndofs+=(v->get())->getNEqs();
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
      for(us i=1;i<geom.nCells-1;i++)
    	vvertex.emplace_back(new TubeVertex());
      // Right *probable* boundary condition
      vvertex.emplace_back(rightTubeVertex());

      us nVertex=vvertex.size();    
      assert(nVertex==geom.nCells);

      for(us i=0;i<vvertex.size();i++){
        vvertex.at(i)->initVertex(i,*this);
      }
      // And initialize again.
      for(us i=0;i<vvertex.size();i++){
        TubeVertex* cvertex=static_cast<TubeVertex*>(vvertex[i].get());
        TRACE(13,"Starting intialization of Vertex "<< i);
        if(i<nVertex-1) cvertex->setRight(*vvertex[i+1].get());
        if(i>0) cvertex->setLeft(*vvertex[i-1].get());
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
      TubeVertex& cvertex=*static_cast<TubeVertex*>(vertex->get());
      mass+=cvertex.rho(0)*cvertex.lg.vVf;
    }
    return mass;
  }

  
  vd Tube::getResAt(us varnr,us freqnr) const{
    TRACE(8,"Tube::getResAt("<<varnr<<","<<freqnr<<")");
    const us& nCells=geom.nCells;
    vd res(nCells);
    assert(varnr<getNDofs());
    
    for(us i=0;i<nCells;i++){
      TubeVertex* cvertex=static_cast<TubeVertex*>(vvertex[i].get());
      switch(varnr) {
        case 0: // Density
          res(i)=cvertex->rho(freqnr);
          break;
        case 1:                 // Volume flown
          res(i)=cvertex->U(freqnr);
          break;
        case 2:                   // Pressure
          res(i)=cvertex->pL()(freqnr);
          break;
        case 3:                 // Temp
          res(i)=cvertex->T()(freqnr);
          break;
        case 4:                 // Temp
          res(i)=cvertex->Ts()(freqnr);
          break;
        default:
          res(i)=0;
      }

    }
    return res;
  }
  vd Tube::getErrorAt(us eqnr,us freqnr) const{
    const us& nCells=geom.nCells;
    vd er(nCells,fillwith::zeros);
    WARN("Unupdated function!!");
    // assert(eqnr<getNDofs());
    // auto eqs=this->getEqs();
    // for(us i=0;i<nCells;i++){
    //   TubeVertex& cvertex=*static_cast<TubeVertex*>(vvertex[i].get());
    //   er(i)=(cvertex.eqs.at(eqnr)->error(cvertex))(freqnr);
    // }
    return er;
  }
  void Tube::resetHarmonics(){
    for(auto v=vvertex.begin();v!=vvertex.end();v++){
      auto &cvertex=*static_cast<TubeVertex*>(v->get());
      cvertex.resetHarmonics();
    }
  }
  void Tube::dmtotdx(vd& dmtotdx_) const{
    TRACE(15,"Tube::dmtotdx()");
    us nvertex=vvertex.size(),Neq;
    us rhodof;
    for(auto v=vvertex.begin();v!=vvertex.end();v++){
      auto &cvertex=*static_cast<TubeVertex*>(v->get());
      rhodof=cvertex.rho.getDofNr();
      dmtotdx_(rhodof)=cvertex.lg.vVf;
    }
  }  
  Tube::~Tube(){
    TRACE(15,"~Tube()");
    cleanup();
  }
  

  
} /* namespace tube */


