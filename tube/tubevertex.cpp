#include "tube.h"
#include "tubevertex.h"
#include "globalconf.h"
#include "weightfactors.h"
#include "jacobian.h"

namespace tube{
  using variable::var;
  using tasystem::Jacobian;

  TubeVertex::TubeVertex(us i,const Tube& tube):
    i(i),
    lg(tube.geom(),i),
    tube(&tube),
    gc(tube.gc),
    c(*this),
    m(*this),
    e(*this),
    sL(*this),
    // s(*this),
    se(*this),
    is(*this)
  {
    vars.resize(5);
    assert(gc!=NULL);

    // Intialize the variables for the right number of harmonics.
    rho_=var(*gc);
    U_=var(*gc);
    T_=var(*gc);
    p_=var(*gc);
    Ts_=var(*gc);

    // Initialize temperature and density variables to something sane
    T_.set(0,gc->T0);
    Ts_.set(0,gc->T0);
    rho_.set(0,gc->rho0());    

    // Fill vars vector
    vars.at(RHONR)=&rho_;
    vars.at(UNR)=&U_;
    vars.at(TNR)=&p_;
    vars.at(PNR)=&T_;
    vars.at(TSNR)=&Ts_;    

    // Fill eqs vector
    // eqs.resize(5);
    eqs.clear();
    eqs.push_back(&c);
    eqs.push_back(&m);
    eqs.push_back(&e);
    // if(i>1 && i<nCells-1)
      // eqs.push_back(&s);
    // else
    eqs.push_back(&sL);
    eqs.push_back(&se);    

  }
  d TubeVertex::getCurrentMass() const{
    return rho_(0)*lg.vVf;
  }
  void TubeVertex::init(const TubeVertex* left,const TubeVertex* right) {
    TRACE(8,"TubeVertex::init(left,right), vertex "<< i << ".");
    // TRACE(25,"Address gc:" <<gc);
    this->left_=left;
    this->right_=right;
    assert(tube);               // *SHOULD* be a valid pointer
    WeightFactors w(*this);
    c.init(w,*tube);
    m.init(w,*tube);
    e.init(w,*tube);
    // s.init(w,*tube);
    sL.init(w,*tube);
    is.init(w,*tube);    

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
  const variable::var& TubeVertex::pR() const {
    TRACE(6,"TubeVertex::pR()");
    assert(right);
    return right_->pL();
  }
  void TubeVertex::setIsentropic(){
    TRACE(15,"TubeVertex::setIsentropic()");
    is.setDofNr(eqs.at(2)->getDofNr());
    eqs[2]=&is;
  }
  
  vd TubeVertex::getp0t() const {
    TRACE(0,"TubeEquation::getp0t()");
    return gc->p0*vd(gc->Ns(),fillwith::ones);
  }    
  vd TubeVertex::errorAt(us eqnr) const{
    TRACE(10,"TubeVertex::errorAt()");
    assert(eqnr<eqs.size());
    return eqs.at(eqnr)->error(*this);
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
      error.subvec(k*Ns,(k+1)*Ns-1)=errorAt(k);
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
  d TubeVertex::getRes(varnr v,us freqnr) const{
    TRACE(4,"TubeVertex::getRes()");
    TRACE(4,"TubeVertex::getRes()");
    switch(v) {
    case varnr::rho: // Density
      return rho()(freqnr);
      break;
    case varnr::U:                 // Volume flown
      return U()(freqnr);
      break;
    case varnr::p:                   // Pressure
      return 0.5*(pL()(freqnr)+pR()(freqnr));
      break;
    case varnr::T:                 // Temp
      return T()(freqnr);
      break;
    case varnr::Ts:                 // Temp
      return Ts()(freqnr);
      break;
    default:
      return 0;
    }
  }
  void TubeVertex::setResVar(varnr v,const vd& res){
    TRACE(4,"TubeVertex::setResVar(varnr,vd)");
    switch(v) {
    case varnr::rho: // Density
      rho_.set(res);
      break;
    case varnr::U:                 // Volume flown
      U_.set(res);
      break;
    case varnr::p:                   // Pressure
      p_.set(res);
      break;
    case varnr::T:                 // Temp
      T_.set(res);
      break;
    case varnr::Ts:                 // Temp
      Ts_.set(res);
      break;
    }
  }
  void TubeVertex::setResVar(varnr v,const variable::var& res){
      TRACE(4,"TubeVertex::setResVar(varnr,var)");
      setResVar(v,res());
  }
  var TubeVertex::getRes(varnr v) const{
    TRACE(4,"TubeVertex::getRes()");
    TRACE(4,"TubeVertex::getRes()");
      switch(v) {
      case varnr::rho: // Density
        return rho();
          break;
      case varnr::U:                 // Volume flown
        return U();
          break;
      case varnr::p:                   // Pressure
          return 0.5*(pL()+pR());
          break;
      case varnr::T:                 // Temp
        return T();
          break;
      case varnr::Ts:                 // Tempc
        return Ts();
          break;
      default:
        return rho();
      }
  }

  void TubeVertex::updateNf(){
    TRACE(10,"TubeVertex::setNf()");
    for(auto var=vars.begin();var!=vars.end();var++)
      (*var)->updateNf();
  }
  void TubeVertex::setRes(const vd& res){
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
  void TubeVertex::show() const{
    cout << "----------------- TubeVertex " << lg.i << "----\n";
    cout << "Showing weight functions for TubeVertex "<< i <<"\n";
    lg.show();
    // w.show();
    // // cout << "cWarti    :"<<cWarti<<"\n";

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

} // namespace tube


