#include "vertex.h"
#include "tube.h"
#include "../var/var.h"
namespace segment{


  Vertex::Vertex(us i,const variable::varoperations& vop1):i(i),vop(vop1),Ns(vop.Ns),rho(vop),U(vop),T(vop),p(vop),Ts(vop){
    TRACE(0,"Vertex constructor");

    vars[0]=&rho;
    vars[1]=&U;
    vars[2]=&T;
    vars[3]=&p;
    vars[4]=&Ts;
   
  }
  Vertex::Vertex(const Vertex& v2):Vertex(v2.i,v2.vop){
    this->rho=v2.rho;
    this->U=v2.U;
    this->T=v2.T;
    this->p=v2.p;
    this->Ts=v2.Ts;
  }
  Vertex&  Vertex::operator=(const Vertex& v2){ // Copy assignment
    this->rho=v2.rho;
    this->U=v2.U;
    this->T=v2.T;
    this->p=v2.p;
    this->Ts=v2.Ts;

    return *this;

  }
  vd Vertex::Error()
  {
    TRACE(0,"Vertex::Error()");
    // TRACE(-1,"Check for position i>0 && i<gp-1...");
    // assert(i>0 && i<tube.geom.gp-1);
    const us Ns=vop.Ns;
    vd error(Neq*Ns);
    for(us k=0;k<Neq;k++)
      {
	error.subvec(k*Ns,(k+1)*Ns-1)=eq[k]->Error();
      }
    return error;
  }
  vd Vertex::GetRes(){			// Get current result vector
    TRACE(0,"Vertex::GetRes()");
    vd res(Neq*Ns);
    for(us k=0;k<Neq;k++){
      res.subvec(k*Ns,k*Ns+Ns-1)=(*vars[k])();
    }
    return res;
  }
  void Vertex::SetRes(vd res){
    TRACE(0,"Vertex::Set()");
    for(us k=0;k<Neq;k++){
      vars[k]->set(res.subvec(k*Ns,k*Ns+Ns-1));
    }
  }
  dmat Vertex::Jac(){		// Return Jacobian
    TRACE(0," Vertex::Jac()...");
    TRACE(-1,"Ns:"<<Ns);
    TRACE(-1,"Neq:"<<Neq);    
    dmat Jac(Neq*Ns,3*Neq*Ns,fillwith::zeros);
    for(us k=0;k<Neq;k++){
      TRACE(-1,"Equation number:"<<k)
	us firstrow=k*Ns;
      // cout << firstrow << " ";
      us firstcol=0;
      // cout << firstcol << " ";
      us lastrow=firstrow+Ns-1;
      // cout << lastrow << " ";
      us lastcol=Jac.n_cols-1;
      // cout << lastcol << " ";
	
      dmat eqJac=eq[k]->Jac();
      // if(k==2)
	// eqJac*=1e6;		// Rescaling of momentum equation
      Jac.submat(firstrow,firstcol,lastrow,lastcol)=eqJac;
      TRACE(-1,"Eqjac returns");

    }
    return Jac;
  }  
  Vertex::~Vertex(){
    TRACE(-5,"Vertex destructor");
  }

} // namespace segment

namespace tube{
  
  TubeVertex::TubeVertex(const Tube& tube1,us i):Vertex(i,tube1.vop),tube(tube1),c(tube,*this),m(tube,*this),e(tube,*this),s(tube,*this),se(tube,*this),is(tube,*this)
  {
    TRACE(0,"TubeVertex contructor");

    eq[0]=&this->c;			// Continuity is first
    eq[1]=&this->m;
    eq[2]=&is; 			// Changed to isentropic
    // eq[2]=&e; 			// Full energy
    eq[3]=&s;
    eq[4]=&se;

    TRACE(0,"TubeVertex constructor done");
  }
  TubeVertex::TubeVertex(const TubeVertex& told):TubeVertex(told.tube,told.i){
    TRACE(0,"TubeVertex::operator(),tgp");
    TRACE(-1,"Copied TubeVertex i:"<<i);
    
  }
  // void TubeVertex::operator=(const TubeVertex& rhs){
  //   TRACE(5,"Error: no copies allowed of TubeVertex");
  //   exit(1);
  // }
  vd  TubeVertex::csource() const {
    TRACE(0,"TubeVertex::csource()");
    return zeros(Ns);}
  vd  TubeVertex::msource() const {
    TRACE(0,"TubeVertex::msource()");
    return zeros(Ns);}
  vd  TubeVertex::esource() const {
    TRACE(0,"TubeVertex::esource()");
    return zeros(Ns);}    
    
    TubeVertex::~TubeVertex(){
    TRACE(-5,"TubeVertex destructor");
  }

 



} // namespace tube
