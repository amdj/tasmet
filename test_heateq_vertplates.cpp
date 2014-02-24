#include "heateq_vertplates.h"
#include <string>

using namespace vertplates;

int main(int argc,char* argv[])
{
	initlog(log4cplus::DEBUG_LOG_LEVEL);
//	cout << log4cplus::FATAL_LOG_LEVEL << endl;

	//initlog((int) *argv[0]);
	LOG4CPLUS_INFO(logger,"Test heateq vertplates started");

	us Nf=5;
	d pm=101325;
	d freq=100;
	d omg=2*pi*freq;
	d T0=293.15;
	d y0=1e-3; //One millimeter

	gas m("air");
	varoperations vop(Nf,freq);
	d rho0=m.rho(T0,pm);
	d cp=m.cp(T0);

	heateq_vertplates h(y0,vop);
//
	var p(vop,pm);
	DEBUGLOG("p:" << p);
	DEBUGLOG("p adata:" << p.getcRes() << endl);
	var T(vop,0.9*T0);
	var Tw(vop,T0);
	var res(vop);
//	cout << "Tw/T:"  << Tw/T;

	var dTdx(vop);
	var dpdx(vop);
	h.setdata(Tw,T,p,dTdx,dpdx,m);

	vd y=linspace(-y0,y0,3);

	vd tempdistr=h.Temp(0,y);
	TRACELOG("y:"<< y);
	DEBUGLOG("Temperature distribution: " << tempdistr);
	cout << tempdistr;
	return 0;
}
