#include "tubevertex.h"
#include "tubeequation.h"
#include "tube.h"
#include "w.h"

namespace tube{


  void TubeVertex::setLeft(const Vertex& v){
    TRACE(8,"TubeVertex::setLeft(vertex)");
    Vertex::setLeft(v);
    assert(this->left==NULL);
    this->left=&static_cast<const TubeVertex&>(v);
  }
  void TubeVertex::setRight(const Vertex& v){
    TRACE(8,"TubeVertex::setRight(vertex)");
    Vertex::setRight(v);
    assert(this->right==NULL);
    this->right=&static_cast<const TubeVertex&>(v);
  }
  void TubeVertex::initTubeVertex(us i,const Tube& thisseg)
  {
    lg=thisseg.geom.localGeom(i);
    TRACE(8,"TubeVertex::initTubeVertex(gc,geom), vertex "<< i << ".");
    // Initialize the Globalconf* ptr and i (the vertex number), 
    Vertex::init(i,*thisseg.gc);	// Which also calls Vertex::updateW()
    // assert(gc!=NULL);
    TRACE(10,"Ns:"<<gc->Ns);
    // Set the zero matrix    
    zero=zeros<dmat>(thisseg.gc->Ns,thisseg.gc->Ns);
    
    // Fill the vector of equation pointers from the Tube instance.
    eq=thisseg.getEq();

    // For compatibility, we store these params in the TubeVertex class.
    nCells=thisseg.geom.nCells;

    // Intialize the variables for the right number of harmonics.
    rho=var(*gc);
    U=var(*gc);
    T=var(*gc);
    p=var(*gc);
    Ts=var(*gc);

    // Initialize density and temperatures
    T.set(0,gc->T0);
    Ts.set(0,gc->T0);
    d T0=gc->T0;
    d p0=gc->p0;
    d rho0=gc->gas.rho(T0,p0);
    rho.set(0,rho0);    
    // Update weight factors
    TRACE(10,"Now running updateW()");
    TubeVertex::updateW(thisseg);
  }
  void TubeVertex::show() const{
    // cout << "----------------- TubeVertex " << lg.i << "----\n";
    cout << "Not up-to-date function.\n";
    cout << "Ts(0):"<< Ts(0) <<"\n";
    w.show();
    cout << "cWddt    :"<<cWddt<<"\n";
    cout << "cWim1    :"<<cWim1<<"\n";
    cout << "cWi      :"<<cWi<<"\n";
    cout << "cWip1    :"<<cWip1<<"\n";
    cout << "cWart    :"<<cWart<<"\n";
    cout << "eWddt    :"<<eWddt<<"\n";
    cout << "eWgim1   :"<<eWgim1<<"\n";
    cout << "eWgi     :"<<eWgi<<"\n";
    cout << "eWgip1   :"<<eWgip1<<"\n";
    cout << "eWkinim1 :"<<eWkinim1<<"\n";
    cout << "eWkini   :"<<eWkini<<"\n";
    cout << "eWkinip1 :"<<eWkinip1<<"\n";
    cout << "eWc1    :"<<eWc1<<"\n";
    cout << "eWc2    :"<<eWc2<<"\n";
    cout << "eWc3    :"<<eWc3<<"\n";
    cout << "eWc4    :"<<eWc4<<"\n";
    
    
  }
  void TubeVertex::updateW(const SegBase& thisseg){
    TRACE(8,"TubeVertex::updateW()");

    const Geom& geom=thisseg.geom;
    const LocalGeom lg=geom.localGeom(i);

    w(*this);			// Weight factors
    updateWEqs(thisseg,w);
  }

  void TubeVertex::updateWEqs(const SegBase& thisseg,const W::W& w){
    TRACE(8,"TubeVertex::updateWEqs()");
    const Geom& geom=thisseg.geom;
    const LocalGeom lg=geom.localGeom(i);
    cWddt=lg.vVf;
    mWddt=lg.vVf/lg.vSf;
    eWddt=lg.vVf;
    eWddtkin=0.5*eWddt/pow(lg.vSf,2);
    cWart=lg.vSf;
    d vSfsq=pow(lg.vSf,2);
    auto& vleft=thisseg.getLeft();
    auto& vright=thisseg.getRight();    
    
    if((i>0 && i<nCells-1) || (i==0 && vleft.size()!=0) || (i==nCells-1 && vright.size()!=0)){
      d vSfLsq=pow(w.vSfL,2);
      d vSfRsq=pow(w.vSfR,2);
      cWim1=-w.UsignL*w.wLl;
      cWi=w.wRl-w.wLr;
      cWip1=w.UsignR*w.wRr;

      mWuim1=-w.UsignL*w.wLl/w.vSfL;
      mWui=w.wRl/lg.vSf-w.wLr/lg.vSf;
      mWuip1=w.UsignR*w.wRr/w.vSfR;

      mWpim1=-w.vSfL*w.wLl;
      mWpi  = lg.vSf*w.wRl-lg.vSf*w.wLr;
      mWpip1= w.vSfR*w.wRr;

      eWgim1=-w.UsignL*w.wLl;
      eWgi=w.wRl-w.wLr;
      eWgip1=w.UsignR*w.wRr;

      eWkinim1=-0.5*w.wLl/vSfLsq;
      eWkini=0.5*w.wRl/vSfsq-w.wLr/vSfLsq;
      eWkinip1=0.5*w.wRr/vSfRsq;
      
      eWc1=-w.vSfL/w.dxm;
      eWc2= lg.vSf/w.dxm;
      eWc3= lg.vSf/w.dxp;
      eWc4=-w.vSfR/w.dxp;
    }
    else if(i==0){
     d vSfRsq=pow(w.vSfR,2);
      cWim1=0;
      cWi=w.wRl;
      cWip1=w.wRr;
      
      mWuim1=0;
      mWui=w.wRl/lg.vSf;
      mWuip1=w.wRr/w.vSfR;
      
      mWpim1=0;
      mWpi=lg.vSf*w.wRl-lg.SfL*w.wL0;
      mWpip1=w.vSfR*w.wRr-lg.SfL*w.wL1;

      eWgim1=0;
      eWgi=w.wRl;
      eWgip1=w.wRr;

      eWkinim1=0;
      eWkini=0.5*w.wRl/vSfsq;
      eWkinip1=0.5*w.wRr/vSfRsq;

      eWc1=0;
      eWc2=0;
      eWc3=lg.vSf/w.dxp;
      eWc4=-w.vSfR/w.dxp;
    }
    else if(i==nCells-1){
      d vSfLsq=pow(w.vSfL,2);
      cWi=-w.wLr;
      cWim1=-w.wLl;
      cWip1=0;

      mWuim1= -w.wLl/lg.SfL;
      mWui=   -w.wLr/lg.SfL;
      mWuip1= 0;
      
      mWpim1=-w.vSfL*w.wLl+lg.SfR*w.wRNm2;
      mWpi=  -lg.vSf*w.wLr+lg.SfR*w.wRNm1;
      mWpip1=0;

      eWgim1=-w.wLl;
      eWgi=-w.wLr;
      eWgip1=0;

      eWkinim1=-0.5*w.wLl/vSfLsq;
      eWkini=-0.5*w.wLr/vSfsq;
      eWkinip1=0;

      eWc1=-w.vSfL/w.dxm;
      eWc2=lg.vSf/w.dxm;
      eWc3=0;
      eWc4=0;
    }
    else{
      WARN("Something went terribly wrong!");
      abort();
    }
    // Contribution from changing cross-sectional area
    mWpi+=lg.SfL-lg.SfR;
  }


    
    
  // void TubeVertex::updateWLeft(const SegBase& thisseg){
  //   int UsignL=1;    
  //   int UsignR=1;    
  //   auto vleft=thisseg.getLeft();
  //   const LocalGeom lg=geom.localGeom(i);
  //   const LocalGeom rlg=geom.localGeom(i+1);

  //   d wRr,wLr;		// Basis weight functions
  //   wRr=wRl=0; 

  //   vxi=lg.vxi;
  //   wLr=(lg.xL-lg.vxim1)/(lg.vxi-lg.vxim1);
  //   vxip1=rlg.vx;
  //   dxp=vxip1-vxi;      
  //   wRr=(lg.xR-lg.vxi)/(lg.vxip1-lg.vxi);
  //   wRl=(lg.vxip1-lg.xR)/(lg.vxip1-lg.vxi);
  //   d vSfLsq=pow(llg.vSf,2);
  //   if(vleft.size()==NULL){
  //     // A left boundary condition can be available, or we have just
  //     // the normal case of an adiabatic wall
  //     TRACE(14,"This TubeVertex is not connected to other vertices on"	\
  // 	    <<" the left side.");
  // 	// We have a boundary on the left side. Weight functions later
  // 	// be customized with derived TubeVertex classes
  //     d vSfsq=pow(lg.vSf,2);
  //     d vSfRsq=pow(rlg.SfR,2);      

  //   }
  //   else{
  //     const SegBase& left=*vleft[0];
  //     TRACE(10,"Jacobian evaluation requires coupling of segments...");
  //     if(left.gettype().compare("Tube")==0){ // Its a Tube
  // 	TRACE(14,"Its a tube!");
  // 	if(left.getRight()[0]->getNumber()==thisseg.getNumber()){
  // 	  TRACE(8,"Connect current head to left segment's tail");
  // 	  const us& LeftnCells=left.geom.nCells;
  // 	  d L=left.geom.L;
  // 	  vxim1=left.geom.vx(LeftnCells-1)-L;
  // 	}
  //     	else{
  // 	  TRACE(8,"Connect current head to left segment's head");
  // 	  TRACE(6,"lg.vxi:"<<lg.vxi);
  // 	  TRACE(6,"lg.vxim1:"<<lg.vxim1);
  // 	  lg.vxim1=-left.geom.vx(0);
  // 	  wLl=(lg.vxi-lg.xL)/(lg.vxi-lg.vxim1);
  // 	  UsignL=-1;
  //     	}
  // 	wLl=(lg.vxi-lg.xL)/(lg.vxi-lg.vxim1);	
  // 	wLr=(lg.xL-lg.vxim1)/(lg.vxi-lg.vxim1);
  //     }
  //     else{			// Its not a Tube
  // 	TRACE(20,"Error, this kind of coupling not yet implemented. Exiting...");
  // 	TRACE(20,"Seg on left side: "<< left.gettype());
  // 	exit(1);
  //     }
  //   } // i==0 and left!=NULL
    



  //   }

    
  // } // updateWLeft()
  // void TubeVertex::updateWRight(const SegBase& thisseg){
  //   int UsignL=1;    
  //   int UsignR=1;    

  //   auto vright=thisseg.getRight();

  //   else if((i==nCells-1) && vright.size()>0) {
  //     TRACE(10,"Jacobian evaluation requires coupling of segments...");
  //     const SegBase& right=*vright[0];
  //     if(right.gettype().compare("Tube")==0){ // Its a Tube
  //   	TRACE(14,"Its a tube!");
  //   	if(right.getLeft()[0]->getNumber()==thisseg.getNumber()){
  //   	  TRACE(8,"Connected current tail to right segment's head");
  //   	  d L=geom.L;
  //   	  lg.vxip1=right.geom.vx(0)+thisseg.geom.L;
  //   	  TRACE(6,"lg.vxi:"<<lg.vxi);
  //   	  TRACE(6,"lg.vxip1:"<<lg.vxip1);
  //   	}
  //     	else{
  //   	  TRACE(8,"Connected current tail to right segment's tail");
  //   	  const us& RightnCells=right.geom.nCells;
  //   	  const d& RightL=right.geom.L;
  //   	  lg.vxip1=RightL-right.geom.vx(nCells)+geom.L;
  //   	  TRACE(6,"lg.vxi:"<<lg.vxi);
  //   	  TRACE(6,"lg.vxip1:"<<lg.vxip1);
  //   	  UsignR=-1;
  //     	}
  //   	wRr=(lg.xR-lg.vxi)/(lg.vxip1-lg.vxi);
  //   	wRl=(lg.vxip1-lg.xR)/(lg.vxip1-lg.vxi);
  //     }
  //     else{			// Its not a Tube
  //   	TRACE(20,"Error, this kind of coupling not yet implemented. Exiting...");
  //   	TRACE(20,"Seg on left side: "<< right.gettype());
  //   	exit(1);
  //     }
  //   } // i==nCells-1 and right!=NULL

  



  //   else if(i==nCells-1 && right!=NULL){
  //     // We have a wall on the right side
  //     TRACE(14,"This TubeVertex is not connected to other vertices on"\
  // 	    <<" the right side.");
  //     LocalGeom llg=geom.localGeom(i-1);
  //     d vSfLsq=pow(llg.SfL,2);
  //     d vSfsq=pow(lg.vSf,2);



  //   } 

  // }
  
  
  vd TubeVertex::getp0t() const {
    TRACE(0,"TubeEquation::getp0t()");
    return gc->p0*vd(gc->Ns,fillwith::ones);
  }    

  vd TubeVertex::error() const
  {
    TRACE(4,"TubeVertex::Error() for TubeVertex "<< i << ".");
    // TRACE(4,"Check for position i>0 && i<gp-1...");
    // assert(i>0 && i<seg.geom.gp-1);
    const us& Ns=gc->Ns;
    TRACE(4,"Assignment of Ns survived:"<< Ns);
    vd error(Neq*Ns);
    for(us k=0;k<Neq;k++)
      {
	error.subvec(k*Ns,(k+1)*Ns-1)=eq[k]->error(*this);
      }
    TRACE(4,"TubeVertex::Error() done.");
    return error;
  }
  vd TubeVertex::getRes() const {			// Get current result vector
    TRACE(4,"TubeVertex::GetRes()");
    const us& Ns=gc->Ns;
    vd res(Neq*Ns);
    for(us k=0;k<Neq;k++){
      res.subvec(k*Ns,k*Ns+Ns-1)=(*vars[k])();
    }
    return res;
  }
  void TubeVertex::setRes(vd res){
    TRACE(4,"TubeVertex::Set()");
    const us& Ns=gc->Ns;
    for(us k=0;k<Neq;k++){
      vars[k]->set(res.subvec(k*gc->Ns,k*Ns+Ns-1));
    }
  }
  dmat TubeVertex::jac() const {		// Return Jacobian
    TRACE(5,"TubeVertex::Jac() for vertex "<< i<< ".");
    const us& Ns=gc->Ns;
    TRACE(5,"Ns:"<<Ns);
    TRACE(5,"Neq:"<<Neq);    
    dmat Jac(Neq*Ns,3*Neq*Ns,fillwith::zeros);
    us firstcol=0;
    us lastcol=Jac.n_cols-1;
    for(us k=0;k<Neq;k++){
      TRACE(5,"Equation number:"<<k);
      us firstrow=k*Ns;
      us lastrow=firstrow+Ns-1;
      Jac.submat(firstrow,firstcol,lastrow,lastcol)=eq[k]->jac(*this);
      TRACE(5,"Equation "<< k <<"... succesfully obtained Jacobian");
    }
    return Jac;
  }  

  
  vd TubeVertex::csource() const {
    TRACE(4,"TubeVertex::csource()");
    return zeros(gc->Ns);}
  vd TubeVertex::msource() const {
    TRACE(4,"TubeVertex::msource()");
    return zeros(gc->Ns);}
  vd TubeVertex::esource() const {
    TRACE(4,"TubeVertex::esource()");
    vd esource=zeros(gc->Ns);
    return esource;
  }    

} // namespace tube
