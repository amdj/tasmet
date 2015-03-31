#include "tube.h"
#include "cell.h"
#include "geom.h"
#include "globalconf.h"
#include "weightfactors.h"
#include "jacobian.h"

namespace tube{
  using variable::var;
  using tasystem::Jacobian;

  Cell::Cell(us i,const Tube& tube):
    tube(&tube),
    gc(&tube.Gc()),
    c(*this),
    m(*this),
    e(*this),
    s(*this),
    // s(*this),
    se(*this),
    is(*this),
    i(i)

  {
    TRACE(15,"Cell::Cell(i,tube)");
    const Geom& geom=tube.geom();
    assert(gc);

    vx=geom.vx(i);
    xR=geom.x(i+1);		// Position of right cell wall
    xL=geom.x(i);			// Position of left cell wall

    // Geometric parameters
    SfL=geom.Sf(i);
    SfR=geom.Sf(i+1);
    SsL=geom.Ss(i);
    SsR=geom.Ss(i+1);

    vSf=geom.vSf(i);
    vSs=geom.vSs(i);
    vVf=geom.vVf(i);
    vVs=geom.vVs(i);
    vrh=geom.vrh(i);

    rhL=geom.rh(i);
    rhR=geom.rh(i+1);
    
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
    vars.resize(constants::nvars);
    vars.at(constants::rho)=&rho_;
    vars.at(constants::rhoU)=&rhoUL_;
    vars.at(constants::T)=&T_;
    vars.at(constants::p)=&p_;
    vars.at(constants::Ts)=&Ts_;    

    // Fill eqs vector
    eqs.resize(constants::neqs);
    eqs.push_back(&c);
    eqs.push_back(&m);
    eqs.push_back(&e);
    eqs.push_back(&s);
    eqs.push_back(&se);    
  }
  Cell::~Cell(){
    TRACE(25,"Cell::~Cell()");
  }
  d Cell::getCurrentMass() const{
    return rho_(0)*vVf;
  }
  void Cell::init(const Cell* left,const Cell* right) {
    TRACE(8,"Cell::init(left,right), cell "<< i << ".");
    // TRACE(25,"Address gc:" <<gc);
    this->left_=left;
    this->right_=right;

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
  d Cell::getValue(Varnr v,us freqnr) const{
    TRACE(4,"Cell::getValue()");
    TRACE(4,"Cell::getValue()");
    switch(v) {
    case Varnr::rho: // Density
      return rho()(freqnr);
      break;
    case Varnr::U:                 // Volume flown
      // return U()(freqnr);
      throw MyError("Not yet implemented");
      break;
    case Varnr::p:                   // Pressure
      return p_()(freqnr);
      break;
    case Varnr::T:                 // Temp
      return T()(freqnr);
      break;
    case Varnr::Ts:                 // Temp
      return Ts()(freqnr);
      break;
    default:
      return 0;
    }
  }
  void Cell::setResVar(Varnr v,const vd& res){
    TRACE(15,"Cell::setResVar(Varnr,vd)");
    switch(v) {
    case Varnr::rho: // Density
      rho_.set(res);
      break;
    case Varnr::U:                 // Volume flown
      WARN("Todo!");
      // U_.set(res);
      break;
    case Varnr::T:                 // Temp
      T_.set(res);
      break;
    case Varnr::p:                   // Pressure
      p_.set(res);
      break;
    case Varnr::Ts:                 // Temp
      Ts_.set(res);
      break;
    default:
      WARN("Varnr" << v << " not handled!");
    }
  }
  void Cell::setResVar(Varnr v,const variable::var& res){
      TRACE(4,"Cell::setResVar(Varnr,var)");
      setResVar(v,res());
  }
  var Cell::getValue(Varnr v) const{
    TRACE(4,"Cell::getValue()");
    TRACE(4,"Cell::getValue()");
      switch(v) {
      case Varnr::rho: // Density
        return rho();
          break;
      case Varnr::rhoU:                 // Volume flown
        WARN("Not good!");
        // return U();
          break;
      case Varnr::p:                   // Pressure
          return p();
          break;
      case Varnr::T:                 // Temp
        return T();
          break;
      case Varnr::Ts:                 // Tempc
        return Ts();
          break;
      }
      // Default:
      return rho();
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
    if(detailnr>=2){
      cout << "Showing weight factors of equations...\n";
      c.show();
      m.show();
      e.show();
      s.show();
    }
    cout <<"Showing LocalGeom data..\n";
    cout <<"i     :" << i<<"\n";
    cout <<"vx    :" << vx<<"\n";
    cout <<"vSf   :" << vSf<<"\n";
    cout <<"SfL   :" << SfL<<"\n";
    cout <<"SfR   :" << SfR<<"\n";    
    cout <<"vVf   :" << vVf<<"\n";
    cout <<"vrh   :" << vrh<<"\n";
    cout <<"xL    :" << xL<<"\n";            
    cout <<"xR    :" << xR<<"\n";
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


