#include "linsolver.h"

namespace solver{

linsolver::linsolver(int gp1,int Nf1): Nf(Nf1),gp(gp1),Ns(2*Nf1+1) {

	//gd=gd::gd(10,3,100.0,"air"); //gp,Nf,freq,mat
	pResult=vd(gp*Ns);
	UResult=vd(gp*Ns);

	DofsPerVar=gp*Ns;
	t=new tube::lintubefv(gd,1,1);
	//cout << "Error size:"  << endl;
	//cout << "Error vector size:" << Error.size() << endl;

}


