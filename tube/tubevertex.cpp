#include "tubevertex.h"
#include "tubeequation.h"
#include "tube.h"
#include "w.h"

namespace tube{

  void TubeVertex::show() const{
    cout << "----------------- TubeVertex " << lg.i << "----\n";
    cout << "Showing weight functions for TubeVertex "<< i <<"\n";
    w.show();
    cout << "cWddt    :"<<cWddt<<"\n";
    cout << "cWim1    :"<<cWim1<<"\n";
    cout << "cWi      :"<<cWi<<"\n";
    cout << "cWip1    :"<<cWip1<<"\n";
    // cout << "cWarti    :"<<cWarti<<"\n";
    cout << "mWuim1   :"<<mWuim1<<"\n";
    cout << "mWui     :"<<mWui<<"\n";
    cout << "mWuip1   :"<<mWuip1<<"\n";

    // cout << "mWpim1   :"<<mWpim1<<"\n";
    // cout << "mWpi     :"<<mWpi<<"\n";
    // cout << "mWpip1   :"<<mWpip1<<"\n";

    cout << "eWddt    :"<<eWddt<<"\n";
    cout << "eWgim1   :"<<eWgim1<<"\n";
    cout << "eWgi     :"<<eWgi<<"\n";
    cout << "eWgip1   :"<<eWgip1<<"\n";
    cout << "eWkinim1 :"<<eWkinim1<<"\n";
    cout << "eWkini   :"<<eWkini<<"\n";
    cout << "eWkinip1 :"<<eWkinip1<<"\n";
    cout << "eWc1     :"<<eWc1<<"\n";
    cout << "eWc2     :"<<eWc2<<"\n";
    cout << "eWc3     :"<<eWc3<<"\n";
    cout << "eWc4     :"<<eWc4<<"\n";

    cout << "Dofnr rho: " << rho.getDofNr() << "\n";
    cout << "Dofnr U  : " << U.getDofNr() << "\n";
    cout << "Dofnr T  : " << T.getDofNr() << "\n";
    cout << "Dofnr p  : " << p.getDofNr() << "\n";
    cout << "Dofnr Ts : " << Ts.getDofNr() << "\n";
    cout << "TubeVertex on left  side:" << left <<"\n";
    cout << "This TubeVertex         :" << this <<"\n";
    cout << "TubeVertex on right side:" << right <<"\n"   ;
  }
  us TubeVertex::getNDofs() const{
    TRACE(5,"TubeVertex::getNDofs()");
    return vars.size()*gc->Ns;
  }
  us TubeVertex::getNEqs() const{
    return eqs.size()*gc->Ns;
  }
  void TubeVertex::setDofNrs(us firstdof){
    TRACE(5,"TubeVertex::setDofNrs()");
    us nvars=vars.size();        // This makes it safe to exclude dofs
                                // in the vars vector
    for(us i=0;i<nvars;i++){
      vars.at(i)->setDofNr(firstdof);
      firstdof+=gc->Ns;
    }
  }
  void TubeVertex::setEqNrs(us firsteq){
    TRACE(5,"TubeVertex::setDofNrs()");
    us neqs=eqs.size();        // This makes it safe to exclude dofs
                                // in the vars vector
    for(us i=0;i<neqs;i++){
      eqs.at(i)->setDofNr(firsteq);
      firsteq+=gc->Ns;
    }
  }  
  void TubeVertex::setLeft(const Vertex& v){
    TRACE(8,"TubeVertex::setLeft(vertex)");
    this->left=&static_cast<const TubeVertex&>(v);
  }
  void TubeVertex::setRight(const Vertex& v){
    TRACE(8,"TubeVertex::setRight(vertex)");
    this->right=&static_cast<const TubeVertex&>(v);
  }
  const variable::var& TubeVertex::pL() const{
    TRACE(6,"TubeVertex::pL()");
    return p;
  }
  const variable::var& TubeVertex::pR() const {
    TRACE(6,"TubeVertex::pR()");
    assert(right);
    return right->p;
  }
  void TubeVertex::initTubeVertex(us i,const Tube& thisseg)
  {
    TRACE(8,"TubeVertex::initTubeVertex(gc,geom), vertex "<< i << ".");
    vars.clear();
    vars.push_back(&rho);
    vars.push_back(&U);
    vars.push_back(&p);
    vars.push_back(&T);
    vars.push_back(&Ts);    
    

    // assert(gc!=NULL);
    TRACE(10,"Ns:"<<gc->Ns);
    // Set the zero matrix    
    zero=zeros<dmat>(thisseg.gc->Ns,thisseg.gc->Ns);
    
    // Fill the vector of equation pointers from the Tube instance.
    auto tubeeqs=thisseg.getEqs();
    eqs.clear(); eqs.reserve(5);
    us eqnr_=0;
    for(auto eq=tubeeqs.begin();eq!=tubeeqs.end();eq++){
      eqs.push_back(std::unique_ptr<TubeEquation>((*eq)->copy()));
      eqs.at(eqnr_)->init(thisseg);
      eqnr_++;
    }

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
    rho.set(0,gc->rho0);    

    // Update weight factors
    TRACE(10,"Now running updateW()");
    TubeVertex::updateW(thisseg);
    // Finally, updating the real weight factors
    updateWEqs(thisseg);
  }
  void TubeVertex::updateW(const SegBase& thisseg){
    TRACE(8,"TubeVertex::updateW()");

    const Geom& geom=thisseg.geom;

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
        connectTubeRight(thisseg);
      }
      else{
        WARN("Right segment's type not understood from connection point of view. Exiting.");
        exit(1);
      }
    }
  }

  void TubeVertex::updateWEqs(const SegBase& thisseg){
    TRACE(8,"TubeVertex::updateWEqs()");
    const Geom& geom=thisseg.geom;
    cWddt=lg.vVf;
    mWddt=lg.vVf/lg.vSf;
    eWddt=lg.vVf;
    eWddtkin=0.5*eWddt/pow(lg.vSf,3);

    d dx=mWddt;
    d vSfsq=pow(w.vSf,2);
    auto& vleft=thisseg.getLeft();
    auto& vright=thisseg.getRight();    
    cWart1=cWart2=cWart3=cWart4=0;		           

    // Always the same
    mWpL=-w.vSf;
    mWpR= w.vSf;

    if(left && right){
      const LocalGeom& llg=left->lg;
      const LocalGeom& rlg=right->lg;      
      d vSfLsq=pow(w.vSfL,2);
      d vSfRsq=pow(w.vSfR,2);
      cWim1=-w.wLl;
      cWi=w.wRl-w.wLr;
      cWip1=w.wRr;

      d vSfLav=0.5*(w.vSf+w.vSfL);
      d vSfRav=0.5*(w.vSf+w.vSfR);

      // This gives numerical issue solver
      // d vSfLav=0.5*(lg.vSf+llg.vSf);
      // d vSfRav=0.5*(lg.vSf+rlg.vSf);
      
      d vSfLavsq=pow(vSfLav,2);
      d vSfRavsq=pow(vSfRav,2);

      // cWart1=-0.5*vSfLav;
      // cWart2= 0.5*vSfLav;
      // cWart3= 0.5*vSfRav;
      // cWart4=-0.5*vSfRav;

      // mWart1=-1;
      // mWart2= 1;
      // mWart3= 1;
      // mWart4=-1;

      // This one should be correct
      mWuim1=-w.wLl/llg.vSf;
      mWui=(w.wRl/lg.vSf-w.wLr/lg.vSf);
      mWuip1=w.wRr/rlg.vSf;

      
      // This gives wiggles but should be correct
      // mWuim1=-w.wLl/vSfLav;
      // mWui=w.wRl/vSfRav-w.wLr/vSfLav;
      // mWuip1=w.wRr/vSfRav;

      // This is *CLEARLY* wrong, but gives nice results. Why???
      // mWuim1=-w.wLl/vSfLav;
      // mWui=w.wRl/vSfLav;
      // mWui+=-w.wLr/vSfRav;
      // mWuip1=w.wRr/vSfRav;

      
      
      // This does not give wiggles, but is wrong! Why???
      // WARN("Wrong but non-wiggly");

      eWgim1=-w.UsignL*w.wLl;
      eWgi=w.wRl-w.wLr;
      eWgip1=w.UsignR*w.wRr;

      eWkinim1=-0.5*w.UsignL*w.wLl/vSfLavsq;
      eWkini=0.5*(w.wRl/vSfRavsq-w.wLr/vSfLavsq);
      eWkinip1=0.5*w.UsignR*w.wRr/vSfRavsq;

      // eWkinim1=-0.5*w.UsignL*w.wLl/vSfLsq;
      // eWkini=0.5*(w.wRl/vSfsq-w.wLr/vSfsq);
      // eWkinip1=0.5*w.UsignR*w.wRr/vSfRsq;
      
      // TRACE(1,"w.dxm:"<< w.dxm);
      // TRACE(1,"w.dxp:"<< w.dxp);
      eWc1=-vSfLav/w.dxm;
      eWc2= vSfLav/w.dxm;
      eWc3= vSfRav/w.dxp;
      eWc4=-vSfRav/w.dxp;

    }
    else if(i==0){
      // Assuming first cell is adiabatic wall
      d vSfRsq=pow(w.vSfR,2);
      cWim1=0;
      cWi=w.wRl;
      cWip1=w.wRr;
      
      mWuim1=0;
      mWui=w.wRl/w.vSf;
      mWuip1=w.wRr/w.vSfR;
      
      eWgim1=0;
      eWgi=w.wRl;
      eWgip1=w.wRr;

      eWkinim1=0;
      eWkini=0.5*w.wRl/vSfsq;
      eWkinip1=0.5*w.wRr/vSfRsq;

      eWc1=0;
      eWc2=0;
      eWc3=w.vSfR/w.dxp;
      eWc4=-w.vSfR/w.dxp;
    }
    else if(i==nCells-1){
      // Assuming last cell is adiabatic wall
      d vSfLsq=pow(w.vSfL,2);
      cWi=-w.wLr;
      cWim1=-w.wLl;
      cWip1=0;

      mWuim1= -w.wLl/w.vSfL;
      mWui=   -w.wLr/w.vSfL;
      mWuip1= 0;
      
      eWgim1=-w.wLl;
      eWgi=-w.wLr;
      eWgip1=0;

      eWkinim1=-0.5*w.wLl/vSfLsq;
      eWkini=-0.5*w.wLr/vSfsq;
      eWkinip1=0;

      eWc1=-w.vSfL/w.dxm;
      eWc2=w.vSfL/w.dxm;
      eWc3=0;
      eWc4=0;
    }
    else{
      WARN("Something went terribly wrong!");
      abort();
    }
    // Contribution from changing cross-sectional area

    // TRACE(50,"Test: put all kin terms 0");
    // eWkinim1=0;
    // eWkini=0;
    // eWkinip1=0;
    // eWddtkin=0;
    // eWgim1=0;
    // eWgi=0;
    // eWgip1=0;
    // eWddt=0;    

  }

  void TubeVertex::connectTubeLeft(const SegBase& thisseg){
    TRACE(15,"TubeVertex::connectTubeLeft()");
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
      TRACE(18,"Forward initializing Tube on left side.");      
      Tube& lefttube_nonconst=const_cast<Tube&>(lefttube);
      lefttube_nonconst.init(*thisseg.gc);
    }

    const us& leftnCells=left.geom.nCells;
    d xvim1;
    if(left.getRight()[0]->getNumber()==thisseg.getNumber()){
      TRACE(8,"Segment " << thisseg.getNumber()<< " connected with "	\
	    << "head to tail of segment"<< left.getNumber() << ".");
      // WE NEED TO BE SURE THAT ALL VERTICES ALREADY HAVE BEEN
      // CREATED. So we initialize the segment from here, if it has
      // not been already.
      this->left=static_cast<const TubeVertex*>(lefttube.vvertex.at(leftnCells-1).get());
      const d& Lleft=left.geom.L;
      xvim1=left.geom.xv(leftnCells-1)-Lleft;
      w.vSfL=left.geom.vSf(leftnCells-1);
    }
    else{
      TRACE(8,"Segment " << thisseg.getNumber()<< " connected with "	\
	    << "head to head of segment"<< left.getNumber() << ".");
      this->left=static_cast<const TubeVertex*>(lefttube.vvertex.at(0).get());
      xvim1=-left.geom.xv(0);
      w.UsignL=-1;
      w.vSfL=left.geom.vSf(0);	
    }
    w.dxm=lg.xvi-xvim1;      
    w.wLl=(lg.xvi)/(lg.xvi-xvim1);	
    w.wLr=1-w.wLl;
  } // connectTubeLeft()
    
    
  void TubeVertex::connectTubeRight(const SegBase& thisseg){
    TRACE(15,"TubeVertex::connectTubeRight()");
    auto vright=thisseg.getRight();
    const Tube& thistube=static_cast<const Tube&>(thisseg);
    const SegBase& right=*thisseg.getRight().at(0);
    const Tube& righttube=static_cast<const Tube&>(right);

    assert(!thistube.bcRight);
    if(righttube.vvertex.size()==0){ // Pre-init the segment
      // For this one, little situation rebuild everything to
      // non-const? I do not think so.
      TRACE(18,"Forward initializing Tube on right side.");
      Tube& righttube_nonconst=const_cast<Tube&>(righttube);
      righttube_nonconst.init(*thisseg.gc);
    }
    
    d xvip1;
    const us& rightnCells=right.geom.nCells;    
    if(right.getLeft()[0]->getNumber()==thisseg.getNumber()){
      TRACE(8,"Connected current tail to right segment's head");
      //   	  d L=geom.L;
      this->right=static_cast<const TubeVertex*>(righttube.vvertex.at(0).get());
      xvip1=right.geom.xv(0)+thisseg.geom.L;
      w.vSfR=right.geom.vSf(0);
    }
    else{
      TRACE(8,"Connected current tail to right segment's tail");
      this->right=static_cast<const TubeVertex*>(righttube.vvertex.at(rightnCells-1).get());
      const d& Lright=right.geom.L;
      xvip1=Lright-right.geom.xv(rightnCells-1)+thisseg.geom.L;
      w.UsignR=-1;
      w.vSfR=right.geom.vSf(rightnCells-1);
    }
    w.dxp=xvip1-lg.xvi;
    w.wRr=(lg.xR-lg.xvi)/(xvip1-lg.xvi);
    w.wRl=(xvip1-lg.xR)/(xvip1-lg.xvi);
  } // connectTubeRight()

  
  
  vd TubeVertex::getp0t() const {
    TRACE(0,"TubeEquation::getp0t()");
    return gc->p0*vd(gc->Ns,fillwith::ones);
  }    

  vd TubeVertex::error() const
  {
    TRACE(4,"TubeVertex::error() for TubeVertex "<< i << ".");
    // TRACE(4,"Check for position i>0 && i<gp-1...");
    // assert(i>0 && i<seg.geom.gp-1);
    const us& Ns=gc->Ns;
    TRACE(4,"Assignment of Ns survived:"<< Ns);
    us Neq=eqs.size();
    us Neqfull=getNEqs();
    vd error(Neqfull);
    for(us k=0;k<Neq;k++){
      error.subvec(k*Ns,(k+1)*Ns-1)=eqs[k]->error(*this);
    }
    TRACE(4,"TubeVertex::error() i="<<i<<" done.");
    return error;
  }
  vd TubeVertex::domg() const
  {
    TRACE(4,"TubeVertex::domg() for TubeVertex "<< i << ".");
    const us& Ns=gc->Ns;
    TRACE(4,"Assignment of Ns survived:"<< Ns);
    us Neq=getNEqs();

    vd domg(Neq*Ns);
    for(us k=0;k<Neq;k++) {
      domg.subvec(k*Ns,(k+1)*Ns-1)=eqs[k]->domg(*this);
    }
    return domg;
  }
  vd TubeVertex::getRes() const {			// Get current result vector
    TRACE(4,"TubeVertex::GetRes()");
    const us& Ns=gc->Ns;
    us nvars=vars.size();        // Only return for number of equations
    vd res(getNDofs());

    for(us k=0;k<nvars;k++){
      res.subvec(k*Ns,k*Ns+Ns-1)=(*vars[k])();
    }
    return res;
  }
  void TubeVertex::setRes(vd res){
    TRACE(10,"TubeVertex::setRes(), i="<< i);
    const us& Ns=gc->Ns;
    us nvars=vars.size();        // Only put in for number of equations
    assert(res.size()==getNDofs());
    for(us k=0;k<nvars;k++){
      vars[k]->set(res.subvec(k*gc->Ns,k*Ns+Ns-1));
    }
    TRACE(10,"TubeVertex::setRes() exiting, i="<< i);    
  }
  void TubeVertex::jac(Jacobian& tofill) const {		// Return Jacobian
    TRACE(5,"TubeVertex::Jac() for vertex "<< i<< ".");
    us neqs=eqs.size();    
    for(us k=0;k<neqs;k++){
      tofill+=eqs[k]->jac(*this);
      TRACE(5,"Equation "<< k <<"... succesfully obtained Jacobian");
    }
    
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
