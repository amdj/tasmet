#include "vertex.h"
#include "var.h"


namespace segment{

  Vertex::Vertex() {
    TRACE(0,"Vertex constructor");
    left=NULL;
    right=NULL;
    vars[0]=&rho;
    vars[1]=&U;
    vars[2]=&T;
    vars[3]=&p;
    vars[4]=&Ts;
   
  }
  Vertex::Vertex(const Vertex& o):Vertex(){}
  Vertex& Vertex::operator=(const Vertex& o) {Vertex(); return *this;}
  void Vertex::Init(us i,const tasystem::Globalconf& g,const Geom& geom){
    TRACE(8,"Vertex::Init()");
    this->i=i;
    this->Ncells=geom.Ncells;
    this->gc=&g;
    this->updateW(geom);
    rho=var(*gc);
    U=var(*gc);
    T=var(*gc);
    p=var(*gc);
    Ts=var(*gc);
  }


  void Vertex::updateW(const Geom& geom)  {
    TRACE(8,"Vertex::updateW()");
    Ncells=geom.Ncells;
    const vd& vx=geom.vx;
    vxi=vx(i);
    vxip1=0;
    vxim1=0;
    // Initialize distances to next node to zero
    dxm=dxp=0;
    if(i>0){
      vxim1=vx(i-1);
      dxm=vxi-vxim1;
    }
    // ****************************** Initalization of vxipm and dxpm
    if(i<Ncells-1){
      vxip1=vx(i+1);
      dxp=vxip1-vxi;
    }
      
    xR=geom.x(i+1);		// Position of right cell wall
    xL=geom.x(i);			// Position of left cell wall
    // Left and right cross-sectional area
    SfL=geom.Sf(i);
    SfR=geom.Sf(i+1);
    // Geometric parameters
    vSf=geom.vSf(i);
    vSs=geom.vSs(i);
    vVf=geom.vVf(i);
    vVs=geom.vVs(i);


    // ****************************** End initialization
    
  }
  // Vertex&  Vertex::operator=(const Vertex& v2){ // Copy assignment
  //   this->rho=v2.rho;
  //   this->U=v2.U;
  //   this->T=v2.T;
  //   this->p=v2.p;
  //   this->Ts=v2.Ts;

  //   return *this;

  // }
  vd Vertex::Error()
  {
    TRACE(7,"Vertex::Error()");
    // TRACE(-1,"Check for position i>0 && i<gp-1...");
    // assert(i>0 && i<seg.geom.gp-1);
    const us& Ns=gc->Ns;
    TRACE(7,"Assignment of Ns survived:"<< Ns);
    vd error(Neq*Ns);
    for(us k=0;k<Neq;k++)
      {
	error.subvec(k*Ns,(k+1)*Ns-1)=eq[k]->Error();
      }
    TRACE(7,"Vertex::Error() done.");
    return error;
  }
  vd Vertex::GetRes(){			// Get current result vector
    TRACE(0,"Vertex::GetRes()");
    const us& Ns=gc->Ns;
    vd res(Neq*Ns);
    for(us k=0;k<Neq;k++){
      res.subvec(k*Ns,k*Ns+Ns-1)=(*vars[k])();
    }
    return res;
  }
  void Vertex::SetRes(vd res){
    TRACE(0,"Vertex::Set()");
    const us& Ns=gc->Ns;
    for(us k=0;k<Neq;k++){
      vars[k]->set(res.subvec(k*gc->Ns,k*Ns+Ns-1));
    }
  }
  dmat Vertex::Jac(){		// Return Jacobian
    const us& Ns=gc->Ns;
    TRACE(0," Vertex::Jac()...");
    TRACE(-1,"Ns:"<<Ns);
    TRACE(-1,"Neq:"<<Neq);    
    dmat Jac(Neq*Ns,3*Neq*Ns,fillwith::zeros);
    for(us k=0;k<Neq;k++){
      TRACE(-1,"Equation number:"<<k);
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

