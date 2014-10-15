#include "geom.h"
#include "grid.h"
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

  void smoothEnds(Geom& first,int firstpos,Geom& second,int secondpos){
    TRACE(3,"SmoothEnd()");

    us i,j;
    
    if(firstpos==FIRST)      {
        i=0;
        TRACE(18,"Smoothing beginning of first Geom, i="<< i);
      }
    else{
      i=first.nCells;
      TRACE(18,"Smoothing end of first Geom, i="<< i);
    }
    if(secondpos==FIRST){
      j=0; 
      TRACE(18,"Smoothing beginning of second Geom, j="<< j);
    }
    else{
      j=second.nCells;
      TRACE(18,"Smoothing end of second Geom, j=" << j);
    }

    vd fx=first.x;
    vd fS=first.S;
    vd fphi=first.phi;
    vd frh=first.rh;
    string fcshape=first.shape;

    vd sx=second.x;
    vd sS=second.S;
    TRACE(1,"Second S:"<< sS);
    vd sphi=second.phi;
    vd srh=second.rh;
    string scshape=second.shape;

    
    // Now adjusting it
    fS(i)=0.5*(first.S(i)+second.S(j));
    sS(j)=0.5*(first.S(i)+second.S(j));

    TRACE(1,"Second S:"<< sS);
    
    fphi(i)=0.5*(first.phi(i)+second.phi(j));
    sphi(j)=0.5*(first.phi(i)+second.phi(j));
    
    first=Geom(fx,fS,fphi,frh,fcshape);
    second=Geom(sx,sS,sphi,srh,scshape);    
  }

  Geom Geom::VertPlates(us gp,d L,d S,d phi,d y0){
    Grid g={gp,L};
    return VertPlates(g,S,phi,y0);
  }
  
  Geom Geom::VertPlates(const Grid& g,d S,d phi,d y0){
    assert(y0>0);
    assert(0<phi && phi<=1.0);
    assert(L>0);

    us gp=g.getgp();
    vd phiv=phi*vd(gp,fillwith::ones);
    vd Sv=S*vd(gp,fillwith::ones);
    vd rh=y0*vd(gp,fillwith::ones);
    return Geom(g,Sv,phiv,rh,"vert");
  }

  
  Geom Geom::Cylinder(us gp,d L,d r){
    return Cone(gp,L,r,r);
  }
  Geom Geom::Cylinder(const Grid& g,d r){
    return Cone(g,r,r);
  }

  Geom Geom::CylinderBlApprox(us gp,d L,d r){
    testgp(gp);
    Grid g={gp,L};
    return CylinderBlApprox(g,r);
  }
  Geom Geom::CylinderBlApprox(const Grid& g,d r){
    TRACE(10,"Geom::CylinderBlApprox()");
    Geom geom=Cone(g,r,r);
    geom.shape="blapprox";
    return geom;
  }  
  Geom Geom::ConeBlApprox(us gp, d L,d r1,d r2){
    testgp(gp);
    Grid g={gp,L};
    return ConeBlApprox(g,r1,r2);
  }
  Geom Geom::ConeBlApprox(const Grid& g,d r1,d r2){
    Geom geom=Cone(g,r1,r2);
    geom.shape="blapprox";
    return geom;
  }

  Geom Geom::Cone(us gp,d L,d r1,d r2){
    TRACE(15,"Geom::Cone()");
    Grid g={gp,L};
    return Cone(g,r1,r2);
  }
  Geom Geom::Cone(const Grid& g,d rleft,d rright){
    TRACE(10,"Geom::Cone()");
    assert(rleft>0);
    assert(rright>0);
    assert(L>0);
    d S1=number_pi*pow(rleft,2);
    d S2=number_pi*pow(rright,2);
    const vd& x1=g.getx();
    const d& L1=g.getL();
    const us gp1=g.getgp();
    TRACE(15,"SFSG");
    vd r1=zeros<vd>(gp1);
    for(us j=0;j<gp1;j++)
      r1(j)=rleft+(rright-rleft)*x1(j)/L1;

    vd phi(gp1,fillwith::ones);
    vd S=number_pi*pow(r1,2);
    vd rh1=r1/2;

    Geom geom(g,S,phi,rh1,"circ");
    if(rleft==rright)
      geom.setPrismatic(true);
    return geom;
  }
  Geom Geom::PrisVertStack(us gp,d L,d S,d phi,d rh){ // Prismatic vertical plates stack
    assert(gp>3);
    assert(L>0);
    assert(S>0);
    assert(phi>0 && phi<1.0);
    Grid g={gp,L};
    return PrisVertStack(g,S,phi,rh);
  }
  Geom Geom::PrisVertStack(const Grid& g,d S,d phi,d rh){ // Prismatic vertical plates stack
    assert(S>0);
    assert(phi>0 && phi<1.0);
    us gp=g.getgp();
    vd phix=phi*vd(gp,fillwith::ones);
    vd Sx=S*vd(gp,fillwith::ones);
    vd rhx=rh*vd(gp,fillwith::ones);
    
    Geom g1(g,Sx,phix,rhx,"vert");
    g1.setPrismatic(true);
    return g1;
  }
  
  d Geom::getFluidVolume() const {
    return arma::sum(vVf);
  }
  Geom::Geom(const Grid& g,const vd& S,const vd& phi,const vd& rh,const string& cshape):Geom(g.getx(),S,phi,rh,cshape){}
  Geom::Geom(const vd& x,const vd& S,const vd& phi,const vd& rh,const string& cshape){
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
    
    xv=vd(nCells);
    vS=vd(nCells);
    vSf=vd(nCells);
    vSs=vd(nCells);
    vrh=vd(nCells);
    vphi=vd(nCells);
    vVf=vd(nCells);
    vVs=vd(nCells);
    for(us j=0;j<nCells;j++){	// Cell-centered scheme
      xv(j)=(x(j+1)+x(j))/2;
      vS(j)=(S(j+1)+S(j))/2;
      vSf(j)=(Sf(j+1)+Sf(j))/2;
      vSs(j)=(Ss(j+1)+Ss(j))/2;
      vrh(j)=(rh(j+1)+rh(j))/2;
      vphi(j)=(phi(j+1)+phi(j))/2;
      vVf(j)=vSf(j)*(x(j+1)-x(j));
      vVs(j)=vSs(j)*(x(j+1)-x(j));
    }
    TRACE(-1,"Celldata xv:"<<xv);
  }
  void Geom::show() const {
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
	     << "rh    : " << rh(0) << "\n"	\
	  ;
      }
    cout << "--------------------------\n"	\
      ;
  }
} // namespace segment
