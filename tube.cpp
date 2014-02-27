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
  // global configuration. Moreover, the tube has gridpoints, "TubeGp" instants. Of these,
  // a tube has gp of them, stored in a vector. In each gridpoint, variables live, which
  // represent the current solution. Moreover, we have equations in each gridpoint. More
  // precisely, in the final solution the continuity, momentum, energy and a suitable
  // equation of state should hold.

namespace segment{
  Seg::Seg(globalconf::Globalconf& g):gc(g){ nL=0; nR=0;} // Seg constructor
}
namespace tube {


  Tube::Tube(globalconf::Globalconf& g,Geom geom):Seg(g),vop(g.Nf,g.freq),geom(geom),gas(gc.gas){
    // Fill vector of gridpoints with data:
    TRACE(5,"Tube constructor started, filling gridpoints vector...");
    for(us i=0; i<geom.gp;i++){
	  TRACE(-1,"i:"<<i);
	  TubeGp ti(*this,i);
      gps.push_back(ti);
    }
    TRACE(5,"Tube constructor done");
    // globalconf instance is put in reference variable gc in inherited class Seg
  }
  Tube::~Tube(){}
  void Tube::Init(d T0,d p0){
    for (std::vector<TubeGp>::iterator it = gps.begin() ; it != gps.end(); ++it){
      it->p.set(p0,0);
      it->T.set(T0,0);
}

}
  Geom::Geom(us gp,d L,d S,d phi,d rh,string cshape):gp(gp){
    x=linspace(0,L,gp);
    this->S=S*ones(gp);
    this->phi=phi*ones(gp);
    this->Sf=phi*this->S;
    this->Ss=(1.0-phi)*this->S;
    this->rh=rh*ones(gp);
    TRACE(-1,"x-vector:" << x);
    TRACE(5,"Simple geom constructor done");
  }
  Geom::~Geom(){}

  TubeGp::TubeGp(const Tube& tube1,us i):tube(tube1),i(i),lc(*this),rho(tube.vop),U(tube.vop),T(tube.vop),p(tube.vop){
    TRACE(0,"TubeGp constructor done");
}
  TubeGp::~TubeGp(){}
  TubeLc::TubeLc(TubeGp& tubegp1): tubegp(tubegp1){
    // TRACE(0,"TubeLc tubegp:"<<tubegp.tube.geom.S);
    TRACE(0,"TubeLc constructor done");
    i=tubegp.i;
    Sf=tubegp.tube.geom.Sf(i);
}
dmat  TubeLc::operator()(){
    // Compute the Jacobian for the subsystem around the current gridpoint
  TRACE(0,"TubeLc::operator()");
    const Geom& geom=tubegp.tube.geom;
    const variable::varoperations& vop=tubegp.tube.vop; // Reference to variable operations
    const Tube& tube=tubegp.tube;
    TRACE(-1,"So far,so good");
    us Ns=vop.Ns;		// Number of samples
    assert(i>0 && i<geom.gp-1);
    TRACE(-1,"So far,so good");
    d dxp=geom.x(i+1)-geom.x(i);
    d dxm=geom.x(i)-geom.x(i-1);
    dmat Uip1td(Ns,Ns,fillwith::zeros);			// Ui+1 in time domain
    dmat Uim1td(Ns,Ns,fillwith::zeros);			// Ui-1 in time domain
    dmat rhoip1td(Ns,Ns,fillwith::zeros);			// rhoi+1 in time domain
    dmat rhoim1td(Ns,Ns,fillwith::zeros);			// rhoi-1 in time domain
    Uip1td.diag()=tube.gps.at(i+1).U.tdata();
    Uim1td.diag()=tube.gps.at(i-1).U.tdata();
    rhoip1td.diag()=tube.gps.at(i+1).rho.tdata();
    rhoim1td.diag()=tube.gps.at(i-1).rho.tdata();

    rhoip1=1.0/(dxp+dxm)*vop.fDFT*Uip1td*vop.iDFT;
    rhoim1=1.0/(dxp+dxm)*vop.fDFT*Uim1td*vop.iDFT;
    Uip1=1.0/(dxp+dxm)*vop.fDFT*rhoip1td*vop.iDFT;
    Uim1=1.0/(dxp+dxm)*vop.fDFT*rhoim1td*vop.iDFT;


    dmat result(Ns,12*Ns,fillwith::zeros);
    // Order is: rho,U,T,p
    TRACE(-1,"Ns:" << Ns);
    TRACE(-1,"rhoim1 size:"<< rhoim1);
    TRACE(-1,"vop dft size:"<< vop.fDFT.size());
    // submat: first row,first col,last row, last col
    result.submat(0,0,Ns-1,Ns-1)=rhoim1;

    return result;
  }
  TubeLc::~TubeLc(){}
  // SubJacobian::SubJacobian(TubeGp& tubegp1):tubegp(tubegp1){
  //   TRACE(0,"Subjacobian constructor");
  //   SubJac.zeros(tubegp.tube.gc.Ns,tubegp.tube.gc.Ns);
  //   TRACE(0,"Subjacobian constructor done");
  // }
  // SubJacobian::~SubJacobian(){}







} /* namespace tube */
