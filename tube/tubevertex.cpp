#include "tubevertex.h"
#include "tubeequation.h"
#include "tube.h"

namespace tube{

  void TubeVertex::show() const{
    cout << "----------------- TubeVertex " << lg.i << "----\n";
    cout << "Showing weight functions for TubeVertex "<< i <<"\n";
    lg.show();
    // w.show();
    // cout << "cWddt    :"<<cWddt<<"\n";
    // cout << "cWim1    :"<<cWim1<<"\n";
    // cout << "cWi      :"<<cWi<<"\n";
    // cout << "cWip1    :"<<cWip1<<"\n";
    // // cout << "cWarti    :"<<cWarti<<"\n";
    // cout << "mWuim1   :"<<mWuim1<<"\n";
    // cout << "mWui     :"<<mWui<<"\n";
    // cout << "mWuip1   :"<<mWuip1<<"\n";

    // cout << "mWpim1   :"<<mWpim1<<"\n";
    // cout << "mWpi     :"<<mWpi<<"\n";
    // cout << "mWpip1   :"<<mWpip1<<"\n";

    // cout << "eWddt    :"<<eWddt<<"\n";
    // cout << "eWgim1   :"<<eWgim1<<"\n";
    // cout << "eWgim    :"<<eWgim<<"\n";
    // cout << "eWgip    :"<<eWgip<<"\n";
    // cout << "eWgip1   :"<<eWgip1<<"\n";
    // cout << "eWgUip1pL:"<<eWgUip1pL<<"\n";
    // cout << "eWgUim1pR:"<<eWgUim1pR<<"\n";    
    
    
    // cout << "eWkinim1 :"<<eWkinim1<<"\n";
    // cout << "eWkini   :"<<eWkini<<"\n";
    // cout << "eWkinip1 :"<<eWkinip1<<"\n";
    // cout << "eWc1     :"<<eWc1<<"\n";
    // cout << "eWc2     :"<<eWc2<<"\n";
    // cout << "eWc3     :"<<eWc3<<"\n";
    // cout << "eWc4     :"<<eWc4<<"\n";

    // cout << "wLl   :"<<wLl<<"\n";
    // cout << "wLr   :"<<wLr<<"\n";
    // cout << "wRl   :"<<wRl<<"\n";
    // cout << "wRr   :"<<wRr<<"\n";
    // cout << "wL0   :"<<wL0<<"\n";
    // cout << "wL1   :"<<wL1<<"\n";
    // cout << "wRNm1 :"<<wRNm1<<"\n";
    // cout << "wRNm2 :"<<wRNm2<<"\n";
    // cout << "vSfL  :"<<vSfL<<"\n";
    // cout << "vSfR  :"<<vSfR<<"\n";
    // cout << "dxm   :"<<dxm<<"\n";
    // cout << "dxp   :"<<dxp<<"\n";
    // cout << "UsignL:"<<UsignL<<"\n";
    // cout << "UsignR:"<<UsignR<<"\n";

    // for(auto eq=eqs.begin();eq!=eqs.end();eq++)
      // (*eq)->show();
    
    // cout << "Number of eqs :" << getNEqs() << "\n";
    // cout << "Number of dofs:" << getNDofs() << "\n";    
    // cout << "Dofnr rho: " << rho.getDofNr() << "\n";
    // cout << "Dofnr U  : " << U.getDofNr() << "\n";
    // cout << "Dofnr p  : " << p.getDofNr() << "\n";
    // cout << "Dofnr T  : " << T.getDofNr() << "\n";
    // cout << "Dofnr Ts : " << Ts.getDofNr() << "\n";
    // cout << "TubeVertex on left  side:" << left <<"\n";
    // cout << "This TubeVertex         :" << this <<"\n";
    // cout << "TubeVertex on right side:" << right <<"\n"   ;
  }
  us TubeVertex::getNDofs() const{
    TRACE(5,"TubeVertex::getNDofs()");
    return vars.size()*gc->Ns();
  }
  us TubeVertex::getNEqs() const{
    return eqs.size()*gc->Ns();
  }
  void TubeVertex::setDofNrs(us firstdof){
    TRACE(5,"TubeVertex::setDofNrs()");
    us nvars=vars.size();        // This makes it safe to exclude dofs
    // in the vars vector
    for(us i=0;i<nvars;i++){
      vars.at(i)->setDofNr(firstdof);
      firstdof+=gc->Ns();
    }
  }
  void TubeVertex::setEqNrs(us firsteq){
    TRACE(5,"TubeVertex::setDofNrs()");
    us neqs=eqs.size();        // This makes it safe to exclude dofs
    // in the vars vector
    for(us i=0;i<neqs;i++){
      eqs.at(i)->setDofNr(firsteq);
      firsteq+=gc->Ns();
    }
  }
  void TubeVertex::resetHarmonics(){
    for(auto var=vars.begin();var!=vars.end();var++)
      (*var)->resetHarmonics();
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
  void TubeVertex::initTubeVertex(us i,const Tube& thistube)
  {
    TRACE(8,"TubeVertex::initTubeVertex(gc,geom), vertex "<< i << ".");

    vars.clear();               // Might be unnessesary
    vars.push_back(&rho);
    vars.push_back(&U);
    vars.push_back(&p);
    vars.push_back(&T);
    vars.push_back(&Ts);    
    

    // assert(gc!=NULL);
    TRACE(10,"Ns:"<<gc->Ns());

    
    // Fill the vector of equation pointers from the Tube instance.
    eqs.clear(); eqs.reserve(6); // Room for one extra equation (minor
                                 // overhead)
    eqs.push_back(&c);
    eqs.push_back(&m);
    eqs.push_back(&e);
    // if(i>1 && i<nCells-1)
      // eqs.push_back(&s);
    // else
    eqs.push_back(&sL);

    
    eqs.push_back(&se);    
    us eqnr_=0;
    for(auto eq=eqs.begin();eq!=eqs.end();eq++){
      (*eq)->init(thistube);
    }

    // For compatibility, we store these params in the TubeVertex class.
    nCells=thistube.geom().nCells();

    // Intialize the variables for the right number of harmonics.
    // TRACE(25,"Address gc:" <<gc);    
    rho=var(*gc);
    U=var(*gc);
    T=var(*gc);
    p=var(*gc);
    Ts=var(*gc);

    // Initialize density and temperatures
    T.set(0,gc->T0);
    Ts.set(0,gc->T0);
    rho.set(0,gc->rho0());    

    // Update weight factors
    TRACE(10,"Now running updateW()");
    TubeVertex::updateW(thistube);
    // Finally, updating the real weight factors
  }
  void TubeVertex::setIsentropic(){
    TRACE(15,"TubeVertex::setIsentropic()");
    is.setDofNr(eqs.at(2)->getDofNr());
    eqs[2]=&is;
  }
  void TubeVertex::updateW(const Tube& thistube){
    TRACE(8,"TubeVertex::updateW()");

    const Geom& geom=thistube.geom();

    vxi=lg.vxi;
    if(i>0) {   
      const LocalGeom& llg=left->lg;
      vxim1=llg.vxi;
      dxm=vxi-vxim1;
      wLl=(lg.vxi-lg.xL)/(lg.vxi-llg.vxi);
      wLr=(lg.xL-llg.vxi)/(lg.vxi-llg.vxi);
      vSfL=llg.vSf;
    }
    if(i==0){
      const LocalGeom& rlg=right->lg;
      vSfL=lg.SfL;
      wL0=rlg.vxi/(rlg.vxi-lg.vxi);
      wL1=-lg.vxi/(rlg.vxi-lg.vxi);
    }
    
    if(i<nCells-1){
      const LocalGeom& rlg=right->lg;
      vxip1=rlg.vxi;
      dxp=vxip1-vxi;      
      vSfR=rlg.vSf;
      wRr=(lg.xR-lg.vxi)/(rlg.vxi-lg.vxi);
      wRl=(rlg.vxi-lg.xR)/(rlg.vxi-lg.vxi);
    }
    if(i==nCells-1){
      const LocalGeom& llg=left->lg;
      wRNm1=(llg.vxi-lg.xR)/(llg.vxi-lg.vxi);
      wRNm2=(lg.xR-lg.vxi)/(llg.vxi-lg.vxi);
      vSfR=lg.SfR;
    }    

    // If we find other segments, we set the left and right pointers
    // to nonzero
    if(i==0 && thistube.getLeft().size()!=0){
      const SegBase& left=*thistube.getLeft().at(0);
      if(left.getType().compare("Tube")==0){ // Its a Tube
        connectTubeLeft(thistube);
        middleVertex();
      }
      else{
        WARN("Left segment's type not understood from connection point of view. Exiting.");
        exit(1);
      }
    }
    else if(i==0){
      leftVertex();
    }
    else if(i==nCells-1 && thistube.getRight().size()!=0){
      const SegBase& right=*thistube.getRight().at(0);
      if(right.getType().compare("Tube")==0){ // Its a Tube
        connectTubeRight(thistube);
        middleVertex();
      }
      else{
        WARN("Right segment's type not understood from connection point of view. Exiting.");
        exit(1);
      }
    }
    else if(i==nCells-1){
      rightVertex();
    }
    else{
      middleVertex();
    }
  } // updateW

    

  void TubeVertex::allVertex(){
    TRACE(5,"TubeVertex::allVertex()");
    
    c.Wddt=lg.vVf;
    m.Wddt=lg.vVf/lg.vSf;
    e.Wddt=lg.vVf;
    e.Wddtkin=0.5*e.Wddt/pow(lg.vSf,2);

    // Always the same
    m.WpL=-lg.vSf;
    m.WpR= lg.vSf;

    if(left){
      sL.WLi=-wLr;
      sL.WLim1=-wLl;
      sL.WLip1=0;
    }
    else{
      sL.WLi=-wL0;
      sL.WLim1=0;
      sL.WLip1=-wL1;
    }
    
  }
  
  void TubeVertex::leftVertex(){
    TRACE(15,"TubeVertex::leftVertex()");
    assert(right);
    allVertex();
    const Geom& geom=*lg.geom;
    const LocalGeom& rlg=right->lg;      

    d& SfL=lg.SfL;
    d& SfR=lg.SfR;
    d SfLsq=pow(SfL,2);
    d SfRsq=pow(SfR,2);
    
    d vSfsq=pow(lg.vSf,2);
    d vSfRsq=pow(vSfR,2);

    c.Wim1=0;
    c.Wi=wRl;
    c.Wip1=wRr;
      
    m.Wuim1=0;
    m.Wui=wRl/lg.vSf;
    m.Wuip1=wRr/vSfR;
      
    e.Wgim1=0;
    e.Wgim=0;
    e.Wgip=wRl;
    e.Wgip1=wRr;

    e.Wkinim1=0;
    e.Wkini=0.5*wRl/vSfsq;
    e.Wkinip1=0.5*wRr/vSfRsq;

    e.Wc1=0;
    e.Wc2=0;
    e.Wc3=vSfR/dxp;
    e.Wc4=-vSfR/dxp;
  }
  void TubeVertex::rightVertex(){
    TRACE(15,"TubeVertex::rightVertex()");
    
    assert(left);
    allVertex();
    const Geom& geom=*lg.geom;
    const LocalGeom& llg=left->lg;

    d& SfL=lg.SfL;
    d SfLsq=pow(SfL,2);
    d vSfsq=pow(lg.vSf,2);
    
    // Assuming last cell is adiabatic wall
    d vSfLsq=pow(vSfL,2);
    c.Wi=-wLr;
    c.Wim1=-wLl;
    c.Wip1=0;

    m.Wuim1= -wLl/vSfL;
    m.Wui=   -wLr/vSfL;
    m.Wuip1= 0;
      
    e.Wgim1=-wLl;
    e.Wgim=-wLr;
    e.Wgip=0;
    e.Wgip1=0;

    e.Wkinim1=-0.5*wLl/vSfLsq;
    e.Wkini=-0.5*wLr/vSfsq;
    e.Wkinip1=0;

    e.Wc1=-SfL/dxm;
    e.Wc2=SfL/dxm;
    e.Wc3=0;
    e.Wc4=0;
  }
  
  void TubeVertex::middleVertex(){
    TRACE(5,"TubeVertex::middleVertex()");
    
    assert(left && right);
    allVertex();

    d vSfsq=pow(lg.vSf,2);

    const LocalGeom& llg=left->lg;
    const LocalGeom& rlg=right->lg;      

    d& SfL=lg.SfL;
    d& SfR=lg.SfR;
    d SfLsq=pow(SfL,2);
    d SfRsq=pow(SfR,2);

    const d& vSfR=rlg.vSf;
    const d& vSfL=llg.vSf;
    d vSfLsq=pow(vSfL,2);
    d vSfRsq=pow(vSfR,2);

    c.Wim1=-wLl;
    c.Wi=wRl-wLr;
    c.Wip1=wRr;

    // d vSfLav=0.5*(lg.vSf+llg.vSf);
    // d vSfRav=0.5*(lg.vSf+rlg.vSf);


    // c.Wart1=-0.5*vSfLav;
    // c.Wart2= 0.5*vSfLav;
    // c.Wart3= 0.5*vSfRav;
    // c.Wart4=-0.5*vSfRav;

    // m.art1=-1;
    // m.art2= 1;
    // m.art3= 1;
    // m.art4=-1;

    // This one should be correct
    m.Wuim1=-wLl/SfL;
    m.Wui=(wRl/SfR-wLr/SfL);
    m.Wuip1=wRr/SfR;

    // But this is also possible
    // m.Wuim1=-wLl/vSfL;
    // m.Wui=(wRl/lg.vSf-wLr/lg.vSf);
    // m.Wuip1=wRr/vSfR;
      
      
    e.Wgim1=-wLl;
    e.Wgim =-wLr;
    e.Wgip = wRl;
    e.Wgip1= wRr;

    // e.Wkinim1=-0.5*UsignL*wLl/vSfLsq;
    // e.Wkini=0.5*(wRl/vSfsq-wLr/vSfsq);
    // e.Wkinip1=0.5*UsignR*wRr/vSfRsq;

    e.Wkinim1=-0.5*UsignL*wLl/SfLsq;
    e.Wkini=0.5*(wRl/SfRsq-wLr/SfLsq);
    e.Wkinip1=0.5*UsignR*wRr/SfRsq;
    

    // TRACE(1,"dxm:"<< dxm);
    // TRACE(1,"dxp:"<< dxp);
    e.Wc1=-SfL/dxm;
    e.Wc2= SfL/dxm;
    e.Wc3= SfR/dxp;
    e.Wc4=-SfR/dxp;

  } 
  
  
  void TubeVertex::connectTubeLeft(const Tube& thistube){
    TRACE(15,"TubeVertex::connectTubeLeft()");
    auto vleft=thistube.getLeft();
    const SegBase& left=*thistube.getLeft().at(0);
    const Tube& lefttube=static_cast<const Tube&>(left);

    // Check for certainty that no internal boundary condition are
    // applied here!
    assert(!thistube.bcLeft);

    if(lefttube.vvertex.size()==0){ // Pre-init the segment
      // For this one, little situation rebuild everything to
      // non-const? I do not think so.
      TRACE(18,"Forward initializing Tube on left side.");      
      Tube& lefttube_nonconst=const_cast<Tube&>(lefttube);
      lefttube_nonconst.init(*thistube.gc);
    }

    const us& leftnCells=left.geom().nCells();
    d vxim1;
    if(left.getRight()[0]->getNumber()==thistube.getNumber()){
      TRACE(8,"Segment " << thistube.getNumber()<< " connected with "	\
            << "head to tail of segment"<< left.getNumber() << ".");
      // WE NEED TO BE SURE THAT ALL VERTICES ALREADY HAVE BEEN
      // CREATED. So we initialize the segment from here, if it has
      // not been already.
      this->left=static_cast<const TubeVertex*>(lefttube.vvertex.at(leftnCells-1).get());
      const d& Lleft=left.geom().L();
      vxim1=left.geom().vx(leftnCells-1)-Lleft;
      vSfL=left.geom().vSf(leftnCells-1);
    }
    else{
      TRACE(8,"Segment " << thistube.getNumber()<< " connected with "	\
            << "head to head of segment"<< left.getNumber() << ".");
      this->left=static_cast<const TubeVertex*>(lefttube.vvertex.at(0).get());
      vxim1=-left.geom().vx(0);
      UsignL=-1;
      vSfL=left.geom().vSf(0);	
    }
    dxm=lg.vxi-vxim1;      
    wLl=(lg.vxi)/(lg.vxi-vxim1);	
    wLr=1-wLl;
  } // connectTubeLeft()
  void TubeVertex::connectTubeRight(const Tube& thistube){
    TRACE(15,"TubeVertex::connectTubeRight()");
    TRACE(15,"SFSG");
    assert(!thistube.getRight().empty());
    auto vright=thistube.getRight();
    const SegBase& right=*vright.at(0);
    TRACE(15,"SFSG");
    const Tube& righttube=static_cast<const Tube&>(right);
    TRACE(15,"SFSG");
    assert(!thistube.bcRight);
    if(righttube.vvertex.size()==0){ // Pre-init the segment
      // For this one, little situation rebuild everything to
      // non-const? I do not think so.
      TRACE(18,"Forward initializing Tube on right side.");
      Tube& righttube_nonconst=const_cast<Tube&>(righttube);
      righttube_nonconst.init(*thistube.gc);
    }
    
    d vxip1;
    const us& rightnCells=right.geom().nCells();    
    if(right.getLeft()[0]->getNumber()==thistube.getNumber()){
      TRACE(8,"Connected current tail to right segment's head");
      //   	  d L=geom().L();
      this->right=static_cast<const TubeVertex*>(righttube.vvertex.at(0).get());
      vxip1=right.geom().vx(0)+thistube.geom().L();
      vSfR=right.geom().vSf(0);
    }
    else{
      TRACE(8,"Connected current tail to right segment's tail");
      this->right=static_cast<const TubeVertex*>(righttube.vvertex.at(rightnCells-1).get());
      const d& Lright=right.geom().L();
      vxip1=Lright-right.geom().vx(rightnCells-1)+thistube.geom().L();
      UsignR=-1;
      vSfR=right.geom().vSf(rightnCells-1);
    }
    dxp=vxip1-lg.vxi;
    wRr=(lg.xR-lg.vxi)/(vxip1-lg.vxi);
    wRl=(vxip1-lg.xR)/(vxip1-lg.vxi);
  } // connectTubeRight()    
    
  
  vd TubeVertex::getp0t() const {
    TRACE(0,"TubeEquation::getp0t()");
    return gc->p0*vd(gc->Ns(),fillwith::ones);
  }    

  vd TubeVertex::error() const
  {
    TRACE(4,"TubeVertex::error() for TubeVertex "<< i << ".");
    // TRACE(4,"Check for position i>0 && i<gp-1...");
    // assert(i>0 && i<seg.geom().gp-1);
    const us& Ns=gc->Ns();
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
  void TubeVertex::domg(vd& domg_) const
  {
    TRACE(4,"TubeVertex::domg() for TubeVertex "<< i << ".");
    const us& Ns=gc->Ns();
    TRACE(4,"Assignment of Ns survived:"<< Ns);
    us neqs=eqs.size();

    for(us k=0;k<neqs;k++) {
      eqs[k]->domg(*this,domg_);
    }
  }
  vd TubeVertex::getRes() const {			// Get current result vector
    TRACE(4,"TubeVertex::GetRes()");
    const us& Ns=gc->Ns();
    us nvars=vars.size();        // Only return for number of equations
    vd res(getNDofs());

    for(us k=0;k<nvars;k++){
      res.subvec(k*Ns,k*Ns+Ns-1)=(*vars[k])();
    }
    return res;
  }
  void TubeVertex::updateNf(){
    TRACE(10,"TubeVertex::setNf()");
    for(auto var=vars.begin();var!=vars.end();var++)
      (*var)->updateNf();
  }
  void TubeVertex::setRes(vd res){
    TRACE(10,"TubeVertex::setRes(), i="<< i);
    const us& Ns=gc->Ns();
    us nvars=vars.size();        // Only put in for number of equations
    assert(res.size()==getNDofs());
    for(us k=0;k<nvars;k++){
      vars[k]->set(res.subvec(k*gc->Ns(),k*Ns+Ns-1));
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
    return zeros(gc->Ns());}
  vd TubeVertex::msource() const {
    TRACE(4,"TubeVertex::msource()");
    return zeros(gc->Ns());}
  vd TubeVertex::esource() const {
    TRACE(4,"TubeVertex::esource()");
    vd esource=zeros(gc->Ns());
    return esource;
  }    

} // namespace tube
