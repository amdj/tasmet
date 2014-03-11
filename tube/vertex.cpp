#include "vertex.h"
#include "tube.h"

namespace tube{

  Vertex::Vertex(us i):i(i){}
  Vertex::~Vertex(){}
  
  TubeVertex::TubeVertex(Tube* tube1,us i):Vertex(i),tube(tube1),Ns(tube->vop.Ns),rho(tube1->vop),U(tube1->vop),T(tube1->vop),
				   p(tube1->vop),Ts(tube1->vop),c(tube1,this),
				   m(tube1,this),e(tube1,this),s(tube1,this),se(tube1,this)
  {
    TRACE(0,"TubeVertex contructor");

    // Geometrical parameters
    vSf=tube->geom.vSf(i);
    vSs=tube->geom.vSs(i);
    vVf=tube->geom.vVf(i);
    vVs=tube->geom.vVs(i);

    wLl=wLr=wRl=wRr=0;		// Initialize weight functions to zero
    d xR=tube->geom.x(i+1);		// Position of right cell wall
    d xL=tube->geom.x(i);		// Position of left cell wall
    vd& vx=tube->geom.vx;
    d& vxi=vx(i);
    const us& Ncells=tube->geom.Ncells;
    d vxip1=0;
    d vxim1=0;
    if(i>0){
      vxim1=vx(i-1);
      wLl=(vxi-xL)/(vxi-vxim1);
      wLr=(xL-vxim1)/(vxi-vxim1);
    }
    if(i<Ncells-1){
      vxip1=vx(i+1);
      wRr=(xR-vxi)/(vxip1-vxi);
      wRl=(vxip1-xR)/(vxip1-vxi);
    }
    

    TRACE(-1,"Testing tube pointer...");
    assert(tube!=NULL);
    eq[0]=&c;			// Continuity is first
    eq[1]=&m;
    eq[2]=&e;
    eq[3]=&s;
    eq[4]=&se;

    vars[0]=&rho;
    vars[1]=&U;
    vars[2]=&T;
    vars[3]=&p;
    vars[4]=&Ts;
    TRACE(0,"TubeVertex constructor done");
  }
  TubeVertex  TubeVertex::operator()(const TubeVertex& tgp){
    TRACE(10,"Error, no copies allowed of a TubeVertex, exiting");
    exit(1);
  }
  vd TubeVertex::Error()
  {
    TRACE(0,"TubeVertex::Error()");
    TRACE(-1,"Check for position i>0 && i<gp-1...");
    assert(i>0 && i<tube->geom.gp-1);
    const us Ns=tube->vop.Ns;
    vd error(Neq*Ns);
    for(us k=0;k<Neq;k++)
      {
	error.subvec(k*Ns,(k+1)*Ns-1)=eq[k]->Error();
      }
    return error;
  }
  vd TubeVertex::Get(){			// Get current result vector
    TRACE(0,"TubeVertex::Get()");
    vd res(Neq*Ns);
    for(us k=0;k<Neq;k++){
      res.subvec(k*Ns,k*Ns+Ns-1)=(*vars[k])();
    }
    return res;
  }
  void TubeVertex::Set(vd res){
    TRACE(0,"TubeVertex::Set()");
    for(us k=0;k<Neq;k++){
      vars[k]->set(res.subvec(k*Ns,k*Ns+Ns-1));
    }
  }
  dmat TubeVertex::Jacobian(){		// Return Jacobian
    TRACE(0,"TubeVertex::Jacobian()");
    TRACE(-1,"Ns:"<<Ns);
    dmat Jac(3*Neq*Ns,3*Neq*Ns);
    for(us k=0;k<Neq;k++){
      Jac.submat(k*Ns,0,k*Ns+Ns-1,3*Neq*Ns-1)=(*eq[k])();
    }
    return Jac;
  }
  void TubeVertex::operator=(const TubeVertex& rhs){
    TRACE(5,"Error: no copies allowed of TubeVertex");
    exit(1);
  }
  TubeVertex::~TubeVertex(){
    TRACE(-5,"TubeVertex destructor");
  }

 



} // namespace tube
