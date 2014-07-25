#include "geom.h"
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
  Geom::Geom(us gp,d L,d Sprismatic,d phiprismatic,d rhprismatic,string shape):shape(shape),gp(gp),Ncells(gp-1){

    TRACE(0,"Geom constructor");
    prismatic=true;
    this->L=L;
    x=linspace(0,L,gp);
 
    S=Sprismatic*ones(gp);
    Ss=(1.0-phiprismatic)*this->S;    
    Sf=phiprismatic*S;

    phi=phiprismatic*ones(gp);
    rh=rhprismatic*ones(gp);

    TRACE(-1,"x-vector:" << x);
    TRACE(5,"Simple geom constructor done");
    Celldata();
  }
  void Geom::Celldata(){
    TRACE(0,"Geom::Celldata");
    vx=vd(Ncells);
    vS=vd(Ncells);
    vSf=vd(Ncells);
    vSs=vd(Ncells);
    vrh=vd(Ncells);
    vphi=vd(Ncells);
    vVf=vd(Ncells);
    vVs=vd(Ncells);
    for(us j=0;j<Ncells;j++){
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
	     << "Ncells: " << Ncells << "\n"	\
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
  Geom::~Geom(){}
} // namespace segment
