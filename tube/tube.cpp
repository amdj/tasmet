/*
 * lintube.cpp

 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#include "tube.h"
#include "tubevertex.h"

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
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
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
    TRACE(14,"Tube::getNDofs()");
    us ndofs=0;
    for(auto v=vvertex.begin();v!=vvertex.end();v++)
      ndofs+=(v->get())->getNDofs();
    return ndofs;
  }
  us Tube::getNEqs() const {
    TRACE(14,"Tube::getNDofs()");
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
    vd rho0=getResAt(0,0);
    return arma::dot(rho0,geom.vVf);
  }

  
  vd Tube::getResAt(us varnr,us freqnr) const{
    TRACE(8,"Tube::getResAt("<<varnr<<","<<freqnr<<")");
    const us& nCells=geom.nCells;
    vd res(nCells);
    assert(varnr<getNDofs());
    for(us i=0;i<nCells;i++){
      TubeVertex* cvertex=static_cast<TubeVertex*>(vvertex[i].get());
      res(i)=cvertex->vars[varnr]->operator()(freqnr);
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
  vd Tube::dmtotdx() const{
    TRACE(15,"Tube::dmtotdx()");
    vd dmtotdx(getNDofs(),fillwith::zeros);
    us nvertex=vvertex.size(),Neq;
    for(us i=0;i<nvertex;i++){
      Neq=vvertex.at(i)->getNDofs();
      dmtotdx(i*Neq*gc->Ns)=geom.vVf(i);
    }
    return dmtotdx;
  }  
  Tube::~Tube(){
    TRACE(15,"~Tube()");
    cleanup();
  }
  

  
} /* namespace tube */


