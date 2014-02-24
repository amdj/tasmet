/* test_lintube.cpp */
#include <complex>
#include <math.h>
#include "lintube.h"
//#include "lintubefv.h"
//#include "freqoperators.h"
//#include "inlinefcn.h"
using namespace std;

int main() {
	us gp=10;
    us Nf=2;
    us Ns=2*Nf+1;
    double f=100;
    double omg=2*pi*f;
    double T=1/f;
    double init=10.0;
    vd sinus(Ns);
    for(us i=0;i<Ns;i++) { sinus(i)=cos(2*pi*i/Ns);}
	vd test=init+sinus;
//
//	var::var pres=var::var(Nf,init);
//	pres.setTdata(test);
//	var::var dpdt=var::var(Nf,init);
//	vd bla=pres.ddt_fd(omg);
//	dpdt.setRes(bla);
//	cout << "d2dt:" << dpdt.ddt_fd(omg) << endl;
	//printVec(pres.getRes());

//	freqoperators::freqoperators fop=freqoperators::freqoperators(Nf,f);
//	cout<< fop.DDTfd << endl;
//	cout<< inv(fop.DDTfd.submat(1,1,Ns-1,Ns-1)) << endl;
	tube::lintube lt=tube::lintube(gp,Nf,1.0,1.0,f);

//	lt.solvesys();
//
	lt.createsys();
	d k=omg/343;
	vd ugok(gp),pgok(gp);
	vd x=lt.getx();
	d L=1.0;
	d up=1.0;
	ugok=up*sin(k*(L-x))/sin(k*L);
//	pyprn("ugok:",ugok);
//	pgok=-1.2*343*up*cos(k*(L-x))/sin(k*L);
//	cout << "pgok:" << pgok << endl;
//	lt.u.setRes(ugok,1);
//	lt.p.setRes(pgok,2);
//	cout<< "p1sin before:" << lt.p.getRes(2) << endl;
//	lt.relax(1000);
	//cout<< "p1cos:" << lt.p.getRes(1) << endl;
//	cout<< "p1sin after:" << lt.p.getRes(2) << endl;
//	cout<< "u1cos after: " << lt.u.getRes(1);
	//cout<< "u1sin: " << lt.u.getRes(2);
	//cout<< lt.fop.iddt*lt.fop.ddt;
//	cout << "Matrix:" << fop.DDTfd << endl;
//	//cout << "Size vector:" << sinus.size() << endl;
//	vd testfd=fop.fDFT*test;
//	vd ddttestFd=fop.DDTfd*fop.DDTfd*testfd;
//	vd ddttestTd=fop.iDFT*ddttestFd;
	//cout<< fop.DDTfd << endl;
	//printVec(ddttestFd);
    return 0;
}
