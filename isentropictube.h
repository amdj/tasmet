#ifndef ISENTROPICTUBE_H
#define ISENTROPICTUBE_H
#include "vtypes.h"
#include "tube.h"


namespace tube {

class isentropictube : public tube {
protected:
	vd Result;
	vd Error;
	us TotalDofs;
    vd x,V;
    double freq,omg;
	const int NumVars=2; //Number of variables involved in problem
	//const double c0=343;
	const double rho0=1.2;
	const double gamm=1.4;
	const double c0sq=gamm*287*293.15;
	double Up;
	const double pm=101325;


public:
    double L,Sf;
    var::var V0,xp,pp;
    isentropictube(us,us,double,double,double,double);
    //			   gp,Nf, L    ,  Sf,  , freq, Up
    virtual ~isentropictube();
	vd getpError(us);
	vd getUError(us);
	vd getpError();
	vd getUError();
	vd getError();
	vd getx();
	void setResult(vd);
	vd getResult();
};

}

#endif // ISENTROPICTUBE_H
