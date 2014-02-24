#include "lintubefv.h"
// Finite volume implementation of linearized 1D Euler equations
namespace tube {

lintubefv::lintubefv(us gp,us Nf,d L,d Sf,d freq,d Up)
 : tube(gp,Nf),L(L),Sf(Sf),freq(freq),Up(Up) {

	x=vd(gp);
	V=vd(gp);
    double h=L/double(gp);
    for(us i=0; i<gp; i++) {
		// Not exactly on the gridpoints, but inbetween
		x(i)=0.5*h+double(i)*h;
		V(i)=h;
    }

    omg=2.0*pi*freq;
	V0=var::var(Nf,h*Sf);
	V0.setRes(Up/omg,2); // Since dVdt=-Up
	TotalDofs=2*DofsPerVar; //For the lintube case, since we only require pressure and velocity.
	Error=vd(TotalDofs);
	Result=vd(TotalDofs);
}


d lintubefv::getError() {
	d e=0;
	return e;
}
vd lintubefv::getx()  {return x; }
lintubefv::~lintubefv() {
    //dtor
}

} // Namespace tube




