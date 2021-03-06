#include "duct.h"
#include "cell.h"
#include "geom.h"
#include "globalconf.h"
#include "jacobian.h"
#include "weightfactors.h"
#include "energy.h"
#include "solidenergy.h"
#include "utils.h"

#define Ns (gc->Ns())

namespace duct{
  using tasystem::var;
  using tasystem::Jacobian;

  Cell::Cell(us i,const Duct& duct):
    duct(&duct),
    gc(&duct.Gc()),
    i(i)
  {
    TRACE(15,"Cell::Cell(i,duct)");
    const Geom& geom=duct.geom();
    assert(gc);

    vx=geom.vx(i);
    xr=geom.x(i+1);		// Position of right cell wall
    xl=geom.x(i);			// Position of left cell wall

    xr=geom.x(i+1);		// Position of right cell wall
    xl=geom.x(i);			// Position of left cell wall

    // Geometric parameters
    Sfl=geom.Sf(i);
    Sfr=geom.Sf(i+1);
    Ssl=geom.Ss(i);
    Ssr=geom.Ss(i+1);

    vSf=geom.vSf(i);
    vSs=geom.vSs(i);
    vVf=geom.vVf(i);
    vVs=geom.vVs(i);
    vrh=geom.vrh(i);

    rhl=geom.rh(i);
    rhr=geom.rh(i+1);
    
    // Intialize the variables for the right number of harmonics.
    rho_=var(*gc);    
    ml_=var(*gc);
    mu_=var(*gc);
    T_=var(*gc);
    p_=var(*gc);
    Tw_=var(*gc);
    Ts_=var(*gc);

    // Initialize temperature and density variables to something sane
    T_.setadata(0,gc->T0());
    Tw_.setadata(0,gc->T0());
    Ts_.setadata(0,gc->T0());
    rho_.setadata(0,gc->rho0());    

  }
  Cell::~Cell(){
    TRACE(15,"Cell::~Cell()");
    utils::purge(eqs);
  }
  d Cell::getMass() const{
    return rho_(0)*vVf;
  }
  void Cell::init(const Cell* left,const Cell* right) {
    TRACE(8,"Cell::init(left,right), cell "<< i << ".");
    // TRACE(25,"Address gc:" <<gc);
    this->left_=left;
    this->right_=right;

    duct->setVarsEqs(*this);
    for (auto& eq : eqs) {
      eq.second->init();
    }
  }
  us Cell::getNDofs() const{
    TRACE(5,"Cell::getNDofs()");
    return vars.size()*Ns;
  }
  us Cell::getNEqs() const{
    return eqs.size()*Ns;
  }
  void Cell::setDofNrs(us firstdof){
    TRACE(5,"Cell::setDofNrs("<<firstdof<<")");
    us nvars=vars.size();        // This makes it safe to exclude dofs
    // in the vars vector
    for (auto& var : vars) { 
      var->setDofNr(firstdof); 
      firstdof+=Ns;
    }
  }
  void Cell::setEqNrs(us firsteq){
    TRACE(5,"Cell::setEqNrs("<<firsteq<<")");
    // in the vars vector
    for (auto& eq : eqs) {
      eq.second->setDofNr(firsteq);
      firsteq+=Ns;
    }
  
  }
  void Cell::resetHarmonics() {
    for(var* v : vars)
      v->resetHarmonics();
  }
  vd Cell::errorAt(us eqnr) const{
    TRACE(10,"Cell::errorAt()");
    // assert(eqnr<eqs.size());
    WARN("Broken function");
    abort();
    // return eqs.at(eqnr)->error();
  }
  vd Cell::error() const
  {
    TRACE(4,"Cell::error() for Cell "<< i << ".");
    // TRACE(4,"Check for position i>0 && i<gp-1...");
    // assert(i>0 && i<seg.geom().gp-1);

    TRACE(4,"Assignment of Ns survived:"<< Ns);
    vd error(getNEqs());
    us k=0;
    for(const auto& eq : eqs){
      error.subvec(k*Ns,(k+1)*Ns-1)=eq.second->error();
      k++;
    }
    TRACE(4,"Cell::error() i="<<i<<" done.");
    return error;
  }
  void Cell::domg(vd& domg_) const
  {
    TRACE(4,"Cell::domg() for Cell "<< i << ".");
    TRACE(4,"Assignment of Ns survived:"<< Ns);
    us neqs=eqs.size();

    for(const auto& eq : eqs)
      eq.second->domg(domg_);
  }
  vd Cell::getRes() const {			// Get current result vector
    TRACE(4,"Cell::GetRes()");
    us nvars=vars.size();        // Only return for number of equations
    vd res(getNDofs());

    for(us k=0;k<nvars;k++){
      // VARTRACE(60,*vars[k]);
      res.subvec(k*Ns,k*Ns+Ns-1)=(*vars[k])();
    }
    return res;
  }
  // vd Cell::U() const {
  //   return (0.5*(rhoUL()+rhoUR())/rho_)();
  // }
  d Cell::getValue(Varnr v,us freqnr) const{
    TRACE(4,"Cell::getValue(Varnr,freqnr)");
    return getValue(v)(freqnr);
  }


  void Cell::setResVar(Varnr v,const vd& res){
    TRACE(15,"Cell::setResVar(Varnr,vd)");
    switch(v) {
    case Varnr::rho: // Density
      rho_.setadata(res);
      break;
    case Varnr::U:                 // Volume flown
      WARN("Todo!");
      // U_.setadata(res);
      break;
    case Varnr::T:                 // Temp
      T_.setadata(res);
      break;
    case Varnr::p:                   // Pressure
      p_.setadata(res);
      break;
    case Varnr::Tw:                 // Temp
      Tw_.setadata(res);
      break;
    case Varnr::Ts:                 // Temp
      Ts_.setadata(res);
      break;
    default:
      WARN("Varnr" << toString(v) << " not handled!");
    }
  }
  void Cell::setResVar(Varnr v,const tasystem::var& res){
      TRACE(4,"Cell::setResVar(Varnr,var)");
      setResVar(v,res());
  }
  var Cell::getValue(Varnr v) const{
    TRACE(4,"Cell::getValue(Varnr)");
      switch(v) {
      case Varnr::none:
        return var(*gc);
      case Varnr::rho: // Density
        return rho();
        break;
      case Varnr::m:                 // Volume flown
        return 0.5*(ml()+mr());
        break;
      case Varnr::mH:                 // Volume flown
        return var(*gc,0.5*(Energy::mHl(*this)+Energy::mHr(*this)));
        break;
      case Varnr::mEkin:
        return var(*gc,0.5*(Energy::mEkinl(*this)+Energy::mEkinr(*this)));	
	break;
      case Varnr::Qs:
      	if(eqs.find(EqType::Sol)!=eqs.end()){
	  const SolidEnergy* e=static_cast<const SolidEnergy*>(eqs.at(EqType::Sol));
          return var(*gc,0.5*(e->QL()+e->QR()));
	}
	else
	  return var(*gc,0);
      case Varnr::F:
	throw MyError("Cell does not contain a force variable.");
      case Varnr::x:
	throw MyError("Cell does not contain a piston position variable.");
      case Varnr::Z:
	throw MyError("Not yet implemented");
      case Varnr::Q:                 // Volume flown
        return var(*gc,0.5*(Energy::QL(*this)+Energy::QR(*this)));
        break;
      case Varnr::p:                   // Pressure
        TRACE(15,"getValue: pressure");
        return p();
        break;
      case Varnr::mu:                   // Pressure
        return mu();
        break;
      case Varnr::T:                 // Temp
        return T();
        break;
      case Varnr::Ts:                 // Temp
        return Ts();
        break;
      case Varnr::Tw:                 // Tempc
        return Tw();
        break;
      case Varnr::U:
        return getValue(Varnr::m)/getValue(Varnr::rho);
        break;
      case Varnr::u:
        return getValue(Varnr::U)/duct->geom().Sf(i);
        break;
      }
      // Default:
      throw MyError("Unknown varnr");
  }
    
  void Cell::updateNf(){
    TRACE(15,"Cell::updateNf()");
    ml_.updateNf();
    rho_.updateNf();
    T_.updateNf();
    p_.updateNf();
    Tw_.updateNf();
    Ts_.updateNf();
    mu_.updateNf();
  }
  void Cell::setRes(const vd& res){
    TRACE(10,"Cell::setRes(), i="<< i);

    us nvars=vars.size();        // Only put in for number of equations
    assert(res.size()==getNDofs());
    us k=0;
    for(var*& v: vars) {
      v->setadata(res.subvec(k*Ns,k*Ns+Ns-1));
      k++;
    }
    TRACE(10,"Cell::setRes() exiting, i="<< i);    
  }
  void Cell::jac(Jacobian& tofill) const {		// Return Jacobian
    TRACE(5,"Cell::Jac() for cell "<< i<< ".");
    us neqs=eqs.size();    
    for(auto &eq : eqs){
      tofill+=eq.second->jac();
      // TRACE(5,"Equation "<< k <<"... succesfully obtained Jacobian");
    }
  }  
  void Cell::show(us detailnr) const{
    cout << "----------------- Cell " << i << "----\n";
    if(detailnr>=4){    
      cout <<"Showing LocalGeom data...\n";
      cout <<"i     :" << i<<"\n";
      cout <<"vx    :" << vx<<"\n";
      cout <<"vSf   :" << vSf<<"\n";
      cout <<"Sfl   :" << Sfl<<"\n";
      cout <<"Sfr   :" << Sfr<<"\n";    
      cout <<"vVf   :" << vVf<<"\n";
      cout <<"vrh   :" << vrh<<"\n";
      cout <<"xl    :" << xl<<"\n";            
      cout <<"xr    :" << xr<<"\n";
    }
    if(detailnr>=5){
      cout << "Showing weight factors of equations...\n";
      WeightFactors(*this).show();
    }
    if(detailnr>=6) {
      cout << "Showing equations...\n";
      for (auto& eq : eqs) {
        eq.second->show();
      }
    } // show

    // cout << "Number of eqs :" << getNEqs() << "\n";
    // cout << "Number of dofs:" << getNDofs() << "\n";    
    // cout << "Dofnr rho: " << rho.getDofNr() << "\n";
    // cout << "Dofnr U  : " << U.getDofNr() << "\n";
    // cout << "Dofnr p  : " << p.getDofNr() << "\n";
    // cout << "Dofnr T  : " << T.getDofNr() << "\n";
    // cout << "Dofnr Tw : " << Tw.getDofNr() << "\n";
    // cout << "Cell on left  side:" << left <<"\n";
    // cout << "This Cell         :" << this <<"\n";
    // cout << "Cell on right side:" << right <<"\n"   ;
  }

} // namespace duct


