/*
 * lintube.cpp

 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#include "tube.h"
using std::cout;
using std::endl;
  // Tried to keep the method definition a bit in order in which a tube is created,
  // including all its components. First a tube is created, which has a geometry and a
  // global configuration. Moreover, the tube has gridpoints, "TubeVertex" instants. Of these,
  // a tube has gp of them, stored in a vector. In each gridpoint, variables live, which
  // represent the current solution. Moreover, we have equations in each gridpoint. More
  // precisely, in the final solution the continuity, momentum, energy and a suitable
  // equation of state should hold.
namespace segment{
  Seg::Seg(globalconf::Globalconf& g):gc(g){ nL=0; nR=0;} // Seg constructor
}
namespace tube {
  Tube::Tube(globalconf::Globalconf& g,Geom geom):Seg(g),vop(g.Nf,g.freq),geom(geom),gas(gc.gas),drag(*this){
    // Fill vector of gridpoints with data:
    TRACE(5,"Tube constructor started, filling gridpoints vector...");
    for(us i=0; i<geom.Ncells;i++){
	  TRACE(-1,"i:"<<i);
	  TubeVertex t(*this,i);
	  vvertex.push_back(t);
    }
    TRACE(5,"Tube constructor done");
    // globalconf instance is put in reference variable gc in inherited class Seg
  }

  vd Tube::Get(){
    TRACE(0,"Tube::Get()");
    const us& Neq=(vvertex[0]).Neq;
    const us& Ns=vop.Ns;
    vd Result(geom.Ncells*vop.Ns*Neq);
    for(us k=0; k<geom.Ncells;k++)
      {
	Result.subvec(k*Ns*Neq,k*Ns*Neq+Ns*Neq-1)=vvertex[k].Get();
      }
    return Result;
  }
  vd Tube::Error(){
    TRACE(0,"Tube::Error(), remember only interior nodes!");
    const us& Neq=(vvertex[0]).Neq;
    const us& Ns=vop.Ns;
    vd error(geom.Ncells*vop.Ns*Neq,fillwith::zeros);
    for(us k=1; k<geom.Ncells-1;k++)
      {
	error.subvec(k*Ns*Neq,k*Ns*Neq+Ns*Neq-1)=vvertex[k].Error();
      }
    return error;
  }
  void Tube::Set(vd res){
    TRACE(0,"Tube::Set");
    const us& Neq=(vvertex[0]).Neq;
    const us& Ns=vop.Ns;
    for(us k=0; k<geom.Ncells;k++)
      {
	vvertex[k].Set(res.subvec(k*Ns*Neq,k*Ns*Neq+Ns*Neq-1));
      }
  }
  void Tube::Init(d T0,d p0){
    for (std::vector<TubeVertex>::iterator it = vvertex.begin() ; it != vvertex.end(); ++it){
      (*it).p.set(p0,0);
      (*it).T.set(T0,0);
      (*it).rho.set(gas.rho(T0,p0),0);
    }
  }
  Tube::~Tube(){
    // for (std::vector<TubeVertex>::iterator it = vvertex.begin() ; it != vvertex.end(); ++it)
    //   {
    //   delete *it;
    //   }
  }
  

} /* namespace tube */

