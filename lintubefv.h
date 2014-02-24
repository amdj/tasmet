#ifndef LINTUBEFV_H
#define LINTUBEFV_H

#include "tube.h"
#include "vvar.h"
#include "vtypes.h"
namespace tube {

class lintubefv : public tube {
protected:
	vd Result;
	vd Error;
	us TotalDofs;
	double Up;

    vd x,V;
	double freq,omg;
	const int NumVars=2; //Number of variables involved in problem
	const double c0=343;
	const double c0sq=pow(c0,2);
	const double pm=101325;
	const double rho0=1.2;
public:
	var::var V0;
    double L,Sf;
    lintubefv(us,us,double,double,double,double);
    virtual ~lintubefv();
	vd getpError(us);
	vd getUError(us);
	vd getpError();
	vd getUError();
	d getError();
	vd getx();
	void setResult(vd);
	vd getResult();


private:
};

}

#endif // LINTUBEFV_H
