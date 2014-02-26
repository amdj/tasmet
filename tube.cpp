/*
 * lintube.cpp

 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#include "tube.h"
using std::cout;
using std::endl;

namespace tube {

tube::tube(us gp,us Nf): gp(gp),Nf(Nf)  {
    Ns=2*Nf+1;
	DofsPerVar=gp*Ns;
	unsigned i;
	//cout << "Nf:" << Nf << endl;
//	rho=vvar::vvar("Density",gp,Nf,1.2);
    p=vvar::vvar("Pressure [Pa]",gp,Nf,100000);
//    U=vvar::vvar("Volume flow",gp,Nf);
    u=vvar::vvar("Velocity",gp,Nf);
//    T=vvar::vvar("Temperature (K)",gp,Nf,293.15);
//    rhoU=vvar::vvar("Mass flow (kg/s)",gp,Nf);
//    rhoUu=vvar::vvar("Momentum flow",gp,Nf);
    //Etot=vvar::vvar("Total energy per unit volume (rho*E)",gp,Nf,1.2*mat->e(293.15));
    //Hf=vvar::vvar("Total enthalpy flow (J/s)",gp,Nf);

}



tube::~tube() {

    unsigned i;

}

} /* namespace tube */
