#include "tube.h"
#include "cell.h"
#include "globalconf.h"
#include "weightfactors.h"
#include "jacobian.h"

namespace tube{
  using variable::var;
  using tasystem::Jacobian;

  Cell::Cell(us i,const Tube& tube):
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
    TRACE(15,"Cell::Cell(i,tube)");

    vars.resize(5);
    assert(gc!=NULL);

    // Intialize the variables for the right number of harmonics.
    rhoUL_=var(*gc);
    rho_=var(*gc);    
    T_=var(*gc);
    p_=var(*gc);
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
  Cell::~Cell(){
    TRACE(25,"Cell::~Cell()");
    delete w_;
  }
  const LocalGeom& Cell::localGeom() const {return weightFactors();}
  const WeightFactors& Cell::weightFactors() const{
    assert(w_);
    return *w_;
  }
  d Cell::getCurrentMass() const{
    return rho_(0)*weightFactors().vVf;
  }
  void Cell::init(const Cell* left,const Cell* right) {
    TRACE(8,"Cell::init(left,right), cell "<< i << ".");
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
  us Cell::getNDofs() const{
    TRACE(5,"Cell::getNDofs()");
    return vars.size()*gc->Ns();
  }
  us Cell::getNEqs() const{
    return eqs.size()*gc->Ns();
  }
  void Cell::setDofNrs(us firstdof){
    TRACE(5,"Cell::setDofNrs()");
    us nvars=vars.size();        // This makes it safe to exclude dofs
    // in the vars vector
    for(us i=0;i<nvars;i++){
      vars.at(i)->setDofNr(firstdof);
      firstdof+=gc->Ns();
    }
  }
  void Cell::setEqNrs(us firsteq){
    TRACE(5,"Cell::setDofNrs()");
    us neqs=eqs.size();        // This makes it safe to exclude dofs
    // in the vars vector
    for(us i=0;i<neqs;i++){
      eqs.at(i)->setDofNr(firsteq);
      firsteq+=gc->Ns();
    }
  }
  void Cell::resetHarmonics() throw(std::exception) {
    for(auto var=vars.begin();var!=vars.end();var++)
      (*var)->resetHarmonics();
  }
  void Cell::setIsentropic(){
    TRACE(15,"Cell::setIsentropic()");
    is.setDofNr(eqs.at(2)->getDofNr());
    eqs[2]=&is;
  }
  
  vd Cell::errorAt(us eqnr) const{
    TRACE(10,"Cell::errorAt()");
    assert(eqnr<eqs.size());
    return eqs.at(eqnr)->error();
  }
  vd Cell::error() const
  {
    TRACE(4,"Cell::error() for Cell "<< i << ".");
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
    TRACE(4,"Cell::error() i="<<i<<" done.");
    return error;
  }
  void Cell::domg(vd& domg_) const
  {
    TRACE(4,"Cell::domg() for Cell "<< i << ".");
    const us& Ns=gc->Ns();
    TRACE(4,"Assignment of Ns survived:"<< Ns);
    us neqs=eqs.size();

    for(us k=0;k<neqs;k++) {
      eqs[k]->domg(domg_);
    }
  }
  vd Cell::getRes() const {			// Get current result vector
    TRACE(4,"Cell::GetRes()");
    const us& Ns=gc->Ns();
    us nvars=vars.size();        // Only return for number of equations
    vd res(getNDofs());

    for(us k=0;k<nvars;k++){
      res.subvec(k*Ns,k*Ns+Ns-1)=(*vars[k])();
    }
    return res;
  }
  vd Cell::U() const {
    return (0.5*(rhoUL()+rhoUR())/rho_)();
  }
  d Cell::getValue(varnr v,us freqnr) const{
    TRACE(4,"Cell::getValue()");
    TRACE(4,"Cell::getValue()");
    switch(v) {
    case varnr::rho: // Density
      return rho()(freqnr);
      break;
    case varnr::U:                 // Volume flown
      // return U()(freqnr);
      throw MyError("Not yet implemented");
      break;
    case varnr::p:                   // Pressure
      return p_;
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
  void Cell::setResVar(varnr v,const vd& res){
    TRACE(15,"Cell::setResVar(varnr,vd)");
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
  void Cell::setResVar(varnr v,const variable::var& res){
      TRACE(4,"Cell::setResVar(varnr,var)");
      setResVar(v,res());
  }
  var Cell::getValue(varnr v) const{
    TRACE(4,"Cell::getValue()");
    TRACE(4,"Cell::getValue()");
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

  void Cell::updateNf(){
    TRACE(10,"Cell::setNf()");
    for(auto var=vars.begin();var!=vars.end();var++)
      (*var)->updateNf();
  }
  void Cell::setRes(const vd& res){
    TRACE(10,"Cell::setRes(), i="<< i);
    const us& Ns=gc->Ns();
    us nvars=vars.size();        // Only put in for number of equations
    assert(res.size()==getNDofs());
    for(us k=0;k<nvars;k++){
      vars[k]->set(res.subvec(k*gc->Ns(),k*Ns+Ns-1));
    }
    TRACE(10,"Cell::setRes() exiting, i="<< i);    
  }
  void Cell::jac(Jacobian& tofill) const {		// Return Jacobian
    TRACE(5,"Cell::Jac() for cell "<< i<< ".");
    us neqs=eqs.size();    
    for(us k=0;k<neqs;k++){
      tofill+=eqs[k]->jac();
      TRACE(5,"Equation "<< k <<"... succesfully obtained Jacobian");
    }
  }  
  vd Cell::csource() const {
    TRACE(4,"Cell::csource()");
    return zeros(gc->Ns());}
  vd Cell::msource() const {
    TRACE(4,"Cell::msource()");
    return zeros(gc->Ns());}
  vd Cell::esource() const {
    TRACE(4,"Cell::esource()");
    vd esource=zeros(gc->Ns());
    return esource;
  }    
  void Cell::show(us detailnr) const{
    cout << "----------------- Cell " << i << "----\n";
    if(detailnr>=4){
      cout << "Showing weight functions for Cell "<< i <<"\n";
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
    // cout << "Cell on left  side:" << left <<"\n";
    // cout << "This Cell         :" << this <<"\n";
    // cout << "Cell on right side:" << right <<"\n"   ;
  }

} // namespace tube


