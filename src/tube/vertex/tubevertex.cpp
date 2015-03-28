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
    tube(&tube),
    gc(&tube.Gc()),
    c(*this),
    m(*this),
    e(*this),
    s(*this),
    // s(*this),
    se(*this),
    is(*this)
  {
    TRACE(15,"TubeVertex::TubeVertex(i,tube)");

    vars.resize(5);
    assert(gc!=NULL);

    // Intialize the variables for the right number of harmonics.
    rho_=var(*gc);
    U_=var(*gc);
    T_=var(*gc);
    pL_=var(*gc);
    Ts_=var(*gc);

    // Initialize temperature and density variables to something sane
    T_.set(0,gc->T0());
    Ts_.set(0,gc->T0());
    rho_.set(0,gc->rho0());    

    // Fill vars vector
    vars.at(RHONR)=&rho_;
    vars.at(UNR)=&U_;
    vars.at(TNR)=&T_;
    vars.at(PNR)=&pL_;
    vars.at(TSNR)=&Ts_;    

    // Fill eqs vector
    // eqs.resize(5);
    eqs.clear();
    eqs.push_back(&c);
    eqs.push_back(&m);
    eqs.push_back(&e);
    eqs.push_back(&s);
    eqs.push_back(&se);    
  }
  TubeVertex::~TubeVertex(){
    delete w_;
  }
  const LocalGeom& TubeVertex::localGeom() const {return weightFactors();}
  const WeightFactors& TubeVertex::weightFactors() const{
    assert(w_);
    return *w_;
  }
  d TubeVertex::getCurrentMass() const{
    return rho_(0)*weightFactors().vVf;
  }
  void TubeVertex::init(const TubeVertex* left,const TubeVertex* right) {
    TRACE(8,"TubeVertex::init(left,right), vertex "<< i << ".");
    // TRACE(25,"Address gc:" <<gc);
    this->left_=left;
    this->right_=right;
    assert(tube);               // *SHOULD* be a valid pointer
    VARTRACE(10,w_);
    delete w_;
    VARTRACE(10,w_);
    w_=new WeightFactors(*this);

    c.init();
    m.init();
    e.init();
    // s.init(w,*tube);
    s.init();
    is.init();    

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
  void TubeVertex::resetHarmonics() throw(std::exception) {
    for(auto var=vars.begin();var!=vars.end();var++)
      (*var)->resetHarmonics();
  }
  void TubeVertex::setIsentropic(){
    TRACE(15,"TubeVertex::setIsentropic()");
    is.setDofNr(eqs.at(2)->getDofNr());
    eqs[2]=&is;
  }
  
  vd TubeVertex::errorAt(us eqnr) const{
    TRACE(10,"TubeVertex::errorAt()");
    assert(eqnr<eqs.size());
    return eqs.at(eqnr)->error();
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
      eqs[k]->domg(domg_);
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
  d TubeVertex::getValue(varnr v,us freqnr) const{
    TRACE(4,"TubeVertex::getValue()");
    TRACE(4,"TubeVertex::getValue()");
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
    TRACE(15,"TubeVertex::setResVar(varnr,vd)");
    switch(v) {
    case varnr::rho: // Density
      rho_.set(res);
      break;
    case varnr::U:                 // Volume flown
      U_.set(res);
      break;
    case varnr::T:                 // Temp
      T_.set(res);
      break;
    case varnr::pL:                   // Pressure
      pL_.set(res);
      break;
    case varnr::Ts:                 // Temp
      Ts_.set(res);
      break;
    default:
      WARN("Varnr" << v << " not handled!");
    }
  }
  void TubeVertex::setResVar(varnr v,const variable::var& res){
      TRACE(4,"TubeVertex::setResVar(varnr,var)");
      setResVar(v,res());
  }
  var TubeVertex::getValue(varnr v) const{
    TRACE(4,"TubeVertex::getValue()");
    TRACE(4,"TubeVertex::getValue()");
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
      tofill+=eqs[k]->jac();
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
  void TubeVertex::show(us detailnr) const{
    cout << "----------------- TubeVertex " << i << "----\n";
    if(detailnr>=4){
      cout << "Showing weight functions for TubeVertex "<< i <<"\n";
      assert(w_);
      w_->show();
    }
    if(detailnr>=2){
      cout << "Showing weight factors of equations...\n";
      c.show();
      m.show();
      e.show();
      s.show();
    }
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


