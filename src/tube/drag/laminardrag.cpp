#include "tube.h"
#include "cell.h"
#include "laminardrag.h"
#include "geom.h"
namespace tube{
  namespace {
    // Resistance force for laminar flow for the zero-frequency.   
    class ZeroFreqDrag{
    public:
      ZeroFreqDrag(const Tube& t);
      d operator()(d mu,d rh,d U) const ; // Full resistance in Newtons
      d operator()(d mu,d rh) const;	// Resistance coefficient
    private:
      d (*zerodrag_funptr)(d,d);
    };

    // Resistance force for laminar flow for the zero-frequency. 
    d zerodrag_vert(d mu,d rh){
      return 3*mu/pow(rh,2);
    }
    d zerodrag_circ(d mu,d rh){
      TRACE(0,"zerodrag_circ");
      return 2*mu/pow(rh,2);
    }
    d zerodrag_blapprox(d mu,d rh){
      TRACE(20,"zerodrag_blapprox");
      return zerodrag_circ(mu,rh);
    }
    d zerodrag_inviscid(d mu,d rh){ return 0; }
    
    ZeroFreqDrag::ZeroFreqDrag(const Tube& t){
      TRACE(0,"ZeroFreqDrag::ZeroFreqDrag()");
      if(t.geom().isBlApprox()){
        TRACE(25,"Blapprox is on");
        zerodrag_funptr=&zerodrag_blapprox;
      }
      else{
        TRACE(25,"Blapprox is off");
        if(t.geom().shape().compare("vert")==0)
          zerodrag_funptr=&zerodrag_vert;
        else if(t.geom().shape().compare("circ")==0){
          TRACE(20,"Circular pore chosen");
          zerodrag_funptr=&zerodrag_circ;
        }
        else if(t.geom().shape().compare("inviscid")==0)
          zerodrag_funptr=&zerodrag_inviscid;
        else
          {
            WARN("Warning: tube.geom.shape unknown for ZeroFreqDrag. Aborting...");
            abort();
          }
      }
    }
    d ZeroFreqDrag::operator()(d mu,d rh,d U) const {
      return (*zerodrag_funptr)(mu,rh)*U;
    }
    d ZeroFreqDrag::operator()(d mu,d rh) const {
      return (*zerodrag_funptr)(mu,rh);
    }

  }// namespace 

  LaminarDragResistance::LaminarDragResistance(const Tube& t)
  {
    TRACE(10,"LaminarDragResistanc::LaminarDragResistance()");
    TRACE(11,"Entering redefinition of Rottfuncs");
    if(t.geom().isBlApprox())
      rf=rottfuncs::RottFuncs("blapprox");
    else
      rf=rottfuncs::RottFuncs(t.geom().shape()); // Reinitialize thermoviscous functions with right shape
    TRACE(11,"Exiting redefinition of Rottfuncs");
    zfd=new ZeroFreqDrag(t);
  }
  LaminarDragResistance::~LaminarDragResistance(){
    delete zfd;
  }
  vd LaminarDragResistance::drag(const Cell& v) const {
    TRACE(10,"LaminarDragResistance::drag(v)");
    vd drag=dm(v)*v.mL()();
    return drag; 		// No momentum scale here, since this is already done in dUi!!!!
  }
  dmat LaminarDragResistance::dm(const Cell& v) const { // Derivative of drag resistance to velocity
    TRACE(10,"LaminarDragResistance::dUi()");
    vc CResistance=ComplexResistancecoef(v);
    variable::var resistance(*v.gc);
    resistance.setadata(CResistance);
    return resistance.freqMultiplyMat();
  }

  vc LaminarDragResistance::ComplexResistancecoef(const Cell& v) const {
    TRACE(0,"LaminarDragResistance::ComplexResistancecoef()");
    const us& Nf=v.gc->Nf();
    const us& i=v.geti();
    d T0=v.T()(0);	// Time-averaged temperature
    d mu0=v.gc->gas().mu(T0);
    d p0=v.p()(0)+v.gc->p0();
    d rho0=v.gc->gas().rho(T0,p0);

    const d& rh=v.vrh;
    d omg=v.gc->getomg();
    vc rescoef(Nf+1);
    rescoef(0)=(*zfd)(mu0,rh);	// Zero frequency drag divided by zero-frequency velocity
    if(Nf>0){
      vd omgvec=omg*linspace(1,Nf,Nf);
      // TRACE(100,"omgvec:\n"<<omgvec);
      vd deltanu=sqrt(2*mu0/(rho0*omgvec));
      vd rh_over_deltanu=rh/deltanu;
      vc fnu=rf.fx(rh/deltanu); // Viscous rott function
      rescoef.subvec(1,Nf)=I*omgvec%(fnu/(1.0-fnu));
    }
    return rescoef;
  }

} // namespace tube








