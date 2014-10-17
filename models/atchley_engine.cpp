#define SMOOTHPERC (1)


#include "models.h"
#include "solver.h"
#include "gas.h"
#include "bc.h"
#include "isentropictube.h"
#include "hopkinslaminarduct.h"
#include "enginesystem.h"
#include "geomhelpers.h"
#include "grid.h"
using namespace segment;
using namespace tasystem;
using namespace tube;
using namespace variable;
using namespace gases;

inline us max(us x,us y){ return x>y?x:y;}

class  Gp{
  d gpfac,dx;
public:
  Gp(d gpfac,d dx):gpfac(gpfac),dx(dx){}
  us operator()(d L){ return max(round(gpfac*L/dx),4);}
};


Solver* Atchley_Engine(d gpfac,us Nf,d freq,d Tr,int loglevel,d kappa,vd p1,d p0,int options)
{
  clearConsole();
  // #Some global params
  inittrace(loglevel);

  // Baseline: 500 gridpoints spread over complete geometry
  // gp: grid refinement factor
  d Ltot=1.0;
  d dx=Ltot/500;
  Gp gp(gpfac,dx);
  
  // d p0=376e3;
  d T0=293.15;  
  Globalconf gc(Nf,freq,"helium",T0,p0,kappa);

  d Rtube=38.2e-3/2;
  // d R1tube=Rtube*1.2;
  d S0=number_pi*pow(Rtube,2);
  d rhtube=Rtube/2;

  // ############################## Resonator
  d Lresendorig=87.97e-2;
  // d Lres=0.774;
  d Lres=Lresendorig;
  // WARN("Wrong length of resonator");
  // d Lres=87.97e-2*1.5;


  Grid resgrid(gp(Lres),Lres);
  // ############################## Cold HX's

  // ############################## Resonator again
  d dxb=1e-5;
  resgrid.setLeftBl(1e-7,1,55*gpfac);
  resgrid.setRightBl(dxb,1,20*gpfac);  
  Geom resgeom=Geom::CylinderBlApprox(resgrid,Rtube);
  // Geom resgeom=Geom::ConeBlApprox(gpres,Lres,R1tube,Rtube);

  // ############################## Cold HX's
  d Lchx=10.2e-3*2;
  d y0chx=1.02e-3/2;
  d phichx=0.70;
  d Lchxgap=1.5e-3;

  Grid chxgrid(gp(Lchx)*10,Lchx);
  chxgrid.setLeftBl(dxb,10,30*gpfac);
  chxgrid.setRightBl(dxb,10,30*gpfac);  
  Geom chxgeom=Geom::VertPlates(chxgrid,S0,phichx,y0chx);

  // ############################## Stack
  d y0stk=0.77e-3/2;
  d Lstk=3.5e-2;
  d phistk=0.73;
  
  Grid stkgrid(gp(Lstk)*10,Lstk);
  stkgrid.setLeftBl(dxb,10,20*gpfac);
  stkgrid.setRightBl(dxb,10,20*gpfac);  
  Geom stkgeom=Geom::VertPlates(stkgrid,S0,phistk,y0stk);

  // ############################## Hot HX
  d y0hhx=y0chx;
  d Lhhx=7.62e-3;
  d phihhx=0.70;
  Geom hhxgeom=Geom::VertPlates(chxgrid,S0,phihhx,y0hhx);

  // ############################## Hot end

  d Lhotend=(5.5e-2)-Lhhx+Lresendorig-Lres;
  Grid hotendgrid(gp(Lhotend)*10,Lhotend);
  hotendgrid.setLeftBl(dxb,25,20*gpfac);
  hotendgrid.setRightBl(1e-7,10,55*gpfac);
  Geom hotendgeom=Geom::CylinderBlApprox(hotendgrid,Rtube);

  
  // And blend togeter
  cout << "Smoothing resonator to chx\n";
  smoothEnds(resgeom,LAST,chxgeom,FIRST,SMOOTHPERC);
  cout << "Smoothing chx to stk\n";
  smoothEnds(chxgeom,LAST,stkgeom,FIRST,SMOOTHPERC);
  cout << "Smoothing hhx to stk\n";  
  smoothEnds(hhxgeom,FIRST,stkgeom,LAST,SMOOTHPERC);
  cout << "Smoothing hotend to hhx\n";  
  smoothEnds(hotendgeom,FIRST,hhxgeom,LAST,SMOOTHPERC);
  
  HopkinsLaminarDuct resonator(resgeom,gc.T0);
  HopkinsLaminarDuct chx(chxgeom,gc.T0);
  HopkinsLaminarDuct stk(stkgeom,gc.T0,Tr);
  HopkinsLaminarDuct hhx(hhxgeom,Tr);
  HopkinsLaminarDuct hotend(hotendgeom,Tr);  


  RightIsoTWall rwall(Tr);
  hotend.addBc(rwall);
  
  if(options & DRIVEN){
    cout << "Driven system\n";
    TaSystem sys(gc);
    var pL(gc);
    for(us i=0;i<gc.Ns();i++)
      pL.set(i,p1(i));
    LeftPressure l1(pL);

    resonator.addBc(l1);
    
    sys.addSeg(resonator);
    sys.addSeg(chx);
    sys.addSeg(stk);
    sys.addSeg(hhx);  
    sys.addSeg(hotend);  
    sys.connectSegs(0,1,tasystem::SegCoupling::tailhead);
    sys.connectSegs(1,2,tasystem::SegCoupling::tailhead);
    sys.connectSegs(2,3,tasystem::SegCoupling::tailhead);
    sys.connectSegs(3,4,tasystem::SegCoupling::tailhead);  
  
    Solver* Sol=new Solver(sys);
    return Sol;
  }
  else{

    EngineSystem sys(gc);
    cout << "NOT driven system\n";
    LeftIsoTWall lwall(T0);
    resonator.addBc(lwall);

    sys.addSeg(resonator);
    sys.addSeg(chx);
    sys.addSeg(stk);
    sys.addSeg(hhx);  
    sys.addSeg(hotend);  
    sys.connectSegs(0,1,tasystem::SegCoupling::tailhead);
    sys.connectSegs(1,2,tasystem::SegCoupling::tailhead);
    sys.connectSegs(2,3,tasystem::SegCoupling::tailhead);
    sys.connectSegs(3,4,tasystem::SegCoupling::tailhead);  
    
    // CHANGED to constraint on second harmonic, see what happens
    sys.setTimingConstraint(0,0,2,2);
    sys.setAmplitudeDof(0,0,2,1);

    // sys.setTimingConstraint(0,0,2,4);
    // sys.setAmplitudeDof(0,0,2,3);
    
    
    Solver* Sol=new Solver(sys);
    return Sol;
  }
}


