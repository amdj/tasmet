#include "tubevertex.h"
#include "tubeequation.h"
#include "tube.h"
#include "w.h"

namespace tube{

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

  void TubeVertex::setLeft(const Vertex& v){
    TRACE(8,"TubeVertex::setLeft(vertex)");
    this->left=&static_cast<const TubeVertex&>(v);
  }
  void TubeVertex::setRight(const Vertex& v){
    TRACE(8,"TubeVertex::setRight(vertex)");
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
  void TubeVertex::updateW(const SegBase& thisseg){
    TRACE(8,"TubeVertex::updateW()");

    const Geom& geom=thisseg.geom;
    const LocalGeom lg=geom.localGeom(i);

    w(*this);			// Weight factors

    // If we find other segments, we set the left and right pointers
    // to nonzero
    if(i==0 && thisseg.getLeft().size()!=0){
      const SegBase& left=*thisseg.getLeft().at(0);
      if(left.getType().compare("Tube")==0){ // Its a Tube
	connectTubeLeft(thisseg);
      }
      else{
	WARN("Left segment's type not understood from connection point of view. Exiting.");
	exit(1);
      }
    }
    if(i==nCells-1 && thisseg.getRight().size()!=0){
      const SegBase& right=*thisseg.getRight().at(0);
      if(right.getType().compare("Tube")==0){ // Its a Tube
	TRACE(100,"connectTubeRight()");
	// connectRight(w);
      }
      else{
	WARN("Right segment's type not understood from connection point of view. Exiting.");
	exit(1);
      }
    }
    // Finally, updating the real weight factors
    updateWEqs(thisseg);
  }

  void TubeVertex::updateWEqs(const SegBase& thisseg){
    TRACE(8,"TubeVertex::updateWEqs()");
    const Geom& geom=thisseg.geom;
    const LocalGeom lg=geom.localGeom(i);
    cWddt=lg.vVf;
    mWddt=lg.vVf/w.vSf;
    eWddt=lg.vVf;
    eWddtkin=0.5*eWddt/pow(w.vSf,2);
    cWart=w.vSf;
    d vSfsq=pow(w.vSf,2);
    auto& vleft=thisseg.getLeft();
    auto& vright=thisseg.getRight();    
    
    if((i>0 && i<nCells-1) || (i==0 && vleft.size()!=0) || (i==nCells-1 && vright.size()!=0)){
      d vSfLsq=pow(w.vSfL,2);
      d vSfRsq=pow(w.vSfR,2);
      cWim1=-w.UsignL*w.wLl;
      cWi=w.wRl-w.wLr;
      cWip1=w.UsignR*w.wRr;

      mWuim1=-w.UsignL*w.wLl/w.vSfL;
      mWui=(w.wRl-w.wLr)/w.vSf;
      mWuip1=w.UsignR*w.wRr/w.vSfR;

      mWpim1=-w.vSfL*w.wLl;
      mWpi  = w.vSf*(w.wRl-w.wLr);
      mWpip1= w.vSfR*w.wRr;

      eWgim1=-w.UsignL*w.wLl;
      eWgi=w.wRl-w.wLr;
      eWgip1=w.UsignR*w.wRr;

      eWkinim1=-0.5*w.wLl/vSfLsq;
      eWkini=0.5*w.wRl/vSfsq-w.wLr/vSfLsq;
      eWkinip1=0.5*w.wRr/vSfRsq;
      
      eWc1=-w.vSfL/w.dxm;
      eWc2= w.vSf/w.dxm;
      eWc3= w.vSf/w.dxp;
      eWc4=-w.vSfR/w.dxp;
    }
    else if(i==0){
      d vSfRsq=pow(w.vSfR,2);
      cWim1=0;
      cWi=w.wRl;
      cWip1=w.wRr;
      
      mWuim1=0;
      mWui=w.wRl/w.vSf;
      mWuip1=w.wRr/w.vSfR;
      
      mWpim1=0;
      mWpi=w.vSf*w.wRl-lg.SfL*w.wL0;
      mWpip1=w.vSfR*w.wRr-lg.SfL*w.wL1;

      eWgim1=0;
      eWgi=w.wRl;
      eWgip1=w.wRr;

      eWkinim1=0;
      eWkini=0.5*w.wRl/vSfsq;
      eWkinip1=0.5*w.wRr/vSfRsq;

      eWc1=0;
      eWc2=0;
      eWc3=w.vSf/w.dxp;
      eWc4=-w.vSfR/w.dxp;
    }
    else if(i==nCells-1){
      d vSfLsq=pow(w.vSfL,2);
      cWi=-w.wLr;
      cWim1=-w.wLl;
      cWip1=0;

      mWuim1= -w.wLl/w.vSfL;
      mWui=   -w.wLr/w.vSfL;
      mWuip1= 0;
      
      mWpim1=-w.vSfL*w.wLl+lg.SfR*w.wRNm2;
      mWpi=  -w.vSf*w.wLr+lg.SfR*w.wRNm1;
      mWpip1=0;

      eWgim1=-w.wLl;
      eWgi=-w.wLr;
      eWgip1=0;

      eWkinim1=-0.5*w.wLl/vSfLsq;
      eWkini=-0.5*w.wLr/vSfsq;
      eWkinip1=0;

      eWc1=-w.vSfL/w.dxm;
      eWc2=w.vSf/w.dxm;
      eWc3=0;
      eWc4=0;
    }
    else{
      WARN("Something went terribly wrong!");
      abort();
    }
    // Contribution from changing cross-sectional area
    mWpi+=w.vSfL-w.vSfR;
  }

  void TubeVertex::connectTubeLeft(const SegBase& thisseg){
    TRACE(15,"TubeVertex::connectLeft()");
    auto vleft=thisseg.getLeft();
    const Tube& thistube=static_cast<const Tube&>(thisseg);
    const SegBase& left=*thisseg.getLeft().at(0);
    const Tube& lefttube=static_cast<const Tube&>(left);

    // Check for certainty that no internal boundary condition are
    // applied here!
    assert(!thistube.bcLeft);

    if(lefttube.vvertex.size()==0){ // Pre-init the segment
      // For this one, little situation rebuild everything to
      // non-const? I do not think so.
      Tube& lefttube_nonconst=const_cast<Tube&>(lefttube);
      lefttube_nonconst.init(*thisseg.gc);
    }

    const us& leftnCells=left.geom.nCells;
    d vxim1;
    if(left.getRight()[0]->getNumber()==thisseg.getNumber()){
      TRACE(8,"Segment " << thisseg.getNumber()<< " connected with "	\
	    << "head to tail of segment"<< left.getNumber() << ".");
      // WE NEED TO BE SURE THAT ALL VERTICES ALREADY HAVE BEEN
      // CREATED. So we initialize the segment from here, if it has
      // not been already.
      this->left=static_cast<const TubeVertex*>(lefttube.vvertex.at(leftnCells-1).get());
      const d& Lleft=left.geom.L;
      vxim1=left.geom.vx(leftnCells-1)-Lleft;
      w.vSfL=left.geom.vSf(leftnCells-1);
    }
    else{
      TRACE(8,"Segment " << thisseg.getNumber()<< " connected with "	\
	    << "head to head of segment"<< left.getNumber() << ".");
      this->left=static_cast<const TubeVertex*>(lefttube.vvertex.at(0).get());
      vxim1=-left.geom.vx(0);
      w.UsignL=-1;
      w.vSfL=left.geom.vSf(0);	
    }
    w.wLl=(lg.vxi)/(lg.vxi-vxim1);	
    w.wLr=1-w.wLl;
  } // connectTubeLeft()
    
    
  void TubeVertex::connectTubeRight(const SegBase& thisseg){
    auto vright=thisseg.getRight();
    const Tube& thistube=static_cast<const Tube&>(thisseg);
    const SegBase& right=*thisseg.getRight().at(0);
    const Tube& righttube=static_cast<const Tube&>(right);

    assert(!thistube.bcRight);
    if(righttube.vvertex.size()==0){ // Pre-init the segment
      // For this one, little situation rebuild everything to
      // non-const? I do not think so.
      Tube& righttube_nonconst=const_cast<Tube&>(righttube);
      righttube_nonconst.init(*thisseg.gc);
    }
    
    d vxip1;
    const us& rightnCells=right.geom.nCells;    
    if(right.getLeft()[0]->getNumber()==thisseg.getNumber()){
      TRACE(8,"Connected current tail to right segment's head");
      //   	  d L=geom.L;
      this->right=static_cast<const TubeVertex*>(righttube.vvertex.at(0).get());
      vxip1=right.geom.vx(0)+thisseg.geom.L;
      w.vSfR=right.geom.vSf(0);
    }
    else{
      TRACE(8,"Connected current tail to right segment's tail");
      this->right=static_cast<const TubeVertex*>(righttube.vvertex.at(rightnCells-1).get());
      const d& Lright=right.geom.L;
      vxip1=Lright-right.geom.vx(rightnCells-1)+thisseg.geom.L;
      w.UsignR=-1;
      w.vSfR=right.geom.vSf(rightnCells-1);
    }
    w.wRr=(lg.xR-lg.vxi)/(vxip1-lg.vxi);
    w.wRl=(vxip1-lg.xR)/(vxip1-lg.vxi);
  } // connectTubeRight()

  
  
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
