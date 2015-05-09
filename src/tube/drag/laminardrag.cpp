#include "tube.h"
#include "cell.h"
#include "laminardrag.h"
#include "geom.h"
#include "math_constants.h"

namespace tube{
  namespace drag {
    using math_common::sq2;
    using rottfuncs::RottFuncs;
  
    // Resistance force for laminar flow for the zero-frequency.   
    inline d mu(const Cell& v) {
      d T0;
      if(v.left()) {
      d wim1,wi; std::tie(wim1,wi)=wf(v);
      T0=wi*v.T()(0)+wim1*v.TL()(0);	// Time-averaged
                                                // temperature at the
                                                // cell wall
      }
      else {
        T0=v.T()(0);	// Time-averaged temperature at the cell wall
      }
      return v.gc->gas().mu(T0);
    }
    inline d rho(const Cell& v) {
      d T0,p0;
      if(v.left()) {
        d wim1,wi; std::tie(wim1,wi)=wf(v);
        T0=wi*v.T()(0)+wim1*v.TL()(0);	// Time-averaged temperature at the cell wall
        p0=v.gc->p0()+wi*v.p()(0)+wim1*v.p()(0);	// Time-averaged temperature at the cell wall
      }
      else {
        T0=v.T()(0);	// Time-averaged temperature at the cell wall
        p0=v.gc->p0()+v.p()(0);	// Time-averaged temperature at the cell wall
      }
      return v.gc->gas().rho(T0,p0);
    }
    // Resistance force for laminar flow for the zero-frequency. 
    d zerodrag_vert(const Cell& v){
      TRACE(5,"zerodrag_vert");
      return 3*mu(v)/(rho(v)*pow(v.rhL,2));
    }
    d zerodrag_circ(const Cell& v){
      TRACE(5,"zerodrag_circ");
      return 2*mu(v)/(rho(v)*pow(v.rhL,2));
    }
    d zerodrag_inviscid(const Cell& v){
      TRACE(5,"zerodrag_inviscid");
      return 0;
    }

    class ZeroFreqDragCoef {
      d (*zerodrag_funptr)(const Cell&);
    public:
      ZeroFreqDragCoef(const Tube& t) {
        TRACE(0,"ZeroFreqDragCoef::ZeroFreqDragCoef()");
        if(t.geom().shape().compare("vert")==0)
          zerodrag_funptr=&zerodrag_vert;
        else if(t.geom().shape().compare("circ")==0){
          TRACE(20,"Circular pore chosen");
          zerodrag_funptr=&zerodrag_circ;
        }
        else if(t.geom().shape().compare("inviscid")==0)
          zerodrag_funptr=&zerodrag_inviscid;
        else {
          WARN("Warning: tube.geom.shape unknown for ZeroFreqDrag. Aborting...");
          abort();
        }
      }
      // We implement a cheap variant of polymorphism
      d operator()(const Cell& v) const {
        return (*zerodrag_funptr)(v);
      }

    };


    LaminarDragResistance::LaminarDragResistance(const Tube& t)
    {
      TRACE(10,"LaminarDragResistanc::LaminarDragResistance()");
      TRACE(11,"Entering redefinition of Rottfuncs");
      if(t.geom().isBlApprox())
        rf=RottFuncs("blapprox");
      else
        rf=RottFuncs(t.geom().shape()); // Reinitialize thermoviscous functions with right shape
      TRACE(11,"Exiting redefinition of Rottfuncs");
      zfd=new ZeroFreqDragCoef(t);
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
    
      const d& rh=v.rhL;

      vc rescoef(Nf+1);
      rescoef(0)=(*zfd)(v);	// Zero frequency drag divided by zero-frequency velocity
      if(Nf>0){
        // Divided by sq2, see Eq. 5.23 of my thesis
        const d omg=v.gc->getomg();
        const vd omgvec=omg*linspace(1,Nf,Nf);
        vd rh_over_deltanu=(shearWaveNumber(v).subvec(1,Nf))/sq2;
        TRACE(20,"SFSF");
        vc fnu=rf.fx(rh_over_deltanu); // Viscous rott function
        TRACE(20,"SFSF");
        rescoef.subvec(1,Nf)=I*(omgvec%(fnu/(1.0-fnu)));
        TRACE(20,"SFSF");
      }
      return rescoef;
    }
  } // namespace drag
} // namespace tube

