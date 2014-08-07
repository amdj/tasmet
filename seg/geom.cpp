#include "geom.h"
#include <assert.h>
namespace segment{
  void testgp(us gp){
    if(gp<4 || gp >MAXGP)
      {
	WARN( "WARNING: Given number of gridpoints is "	\
	     << gp<<", which is larger than MAXGP.\n"
	      << "MAXGP is: " << MAXGP << ". Now exiting...");
	exit(1);
	  
      }
  }
  Geom Geom::Cylinder(us gp, d L,d r){
    return Cone(gp,L,r,r);
  }
  Geom Geom::Cone(us gp,d L,d r1,d r2){
    TRACE(10,"Geom::Cone()");
    assert(gp>3);
    assert(r1>0);
    assert(r2>0);
    assert(L>0);
    d S1=number_pi*pow(r1,2);
    d S2=number_pi*pow(r2,2);
    vd r=linspace(r1,r2,gp);
    vd x=linspace(0,L,gp);
    vd phi(gp,fillwith::ones);
    vd S=number_pi*pow(r,2);
    vd rh=S/(2.0*number_pi*r);

    Geom geom(x,S,phi,rh,"circ");
    if(r1==r2)
      geom.setPrismatic(true);
    return geom;
  }
  Geom Geom::PrisVertStack(us gp,d L,d S,d phi,d rh){ // Prismatic vertical plates stack
    assert(gp>3);
    assert(L>0);
    assert(S>0);
    assert(phi>0 && phi<1.0);

    vd x=linspace(0,L,gp);
    vd phix=phi*vd(gp,fillwith::ones);
    vd Sx=S*vd(gp,fillwith::ones);
    vd rhx=rh*vd(gp,fillwith::ones);
    
    Geom g1(x,Sx,phix,rhx,"vert");
    g1.setPrismatic(true);
    return g1;
  }
  Geom::Geom(vd& x,vd& S,vd& phi,vd& rh,string cshape){
    TRACE(10,"Geom constructor");
    // Sanity checks
    assert(max(phi)<=1.0);
    assert(min(phi)>=0);
    assert(min(S)>0);
    assert(min(x)>=0);
    assert(x.size()>0);
    this->x=x;
    this->shape=cshape;
    this->S=S;
    this->phi=phi;
    this->rh=rh;
    if(min(S)==max(S))
      setPrismatic(true);
    Sf=phi%S;
    Ss=(1.0-phi)%S;

    this->L=x(x.size()-1);
    gp=x.size();
    nCells=gp-1;
    
    vx=vd(nCells);
    vS=vd(nCells);
    vSf=vd(nCells);
    vSs=vd(nCells);
    vrh=vd(nCells);
    vphi=vd(nCells);
    vVf=vd(nCells);
    vVs=vd(nCells);
    for(us j=0;j<nCells;j++){	// Cell-centered scheme
      vx(j)=(x(j+1)+x(j))/2;
      vS(j)=(S(j+1)+S(j))/2;
      vSf(j)=(Sf(j+1)+Sf(j))/2;
      vSs(j)=(Ss(j+1)+Ss(j))/2;
      vrh(j)=(rh(j+1)+rh(j))/2;
      vphi(j)=(phi(j+1)+phi(j))/2;
      vVf(j)=vSf(j)*(x(j+1)-x(j));
      vVs(j)=vSs(j)*(x(j+1)-x(j));
    }
    TRACE(-1,"Celldata vx:"<<vx);
  }
  void Geom::show(){
    if(true)
      {
	cout << "-------- Geometry --------\n"	\
	     << "nCells: " << nCells << "\n"	\
	     << "L     : " << L<< "\n"		\
	     << "Shape : " << shape<< "\n"	\
	  ;
      }
    if(prismatic)
      {
	cout << "S     : " << S(0) << "\n"	\
	     << "Sf    : " << Sf(0) << "\n"	\
	     << "phi   : " << phi(0) << "\n"	\
	  ;
      }
    cout << "--------------------------\n"	\
      ;
  }
} // namespace segment
