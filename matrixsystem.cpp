#include "matrixsystem.h"

matrixsystem::matrixsystem(us nvars,us Nf,us gp): gp(gp),Nf(Nf),nvars(nvars) {
	Ns=2*Nf+1;
	TotalDofs=nvars*Ns*gp;
	Kglob=zeros<dmat>(TotalDofs,TotalDofs);
	Rglob=zeros<vd>(TotalDofs);
}

matrixsystem::setsubsys(us eqnr,us gp,us varj,us gpk){ //Matrix dLi_dpk: derivative
// of equation i to variable varj on gridpoint k.


}

matrixsystem::~matrixsystem()
{
	//dtor
}
