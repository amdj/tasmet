#include "linsolver.h"
#include "lintube.h"

using namespace solver;


int main() {
	linsolver a = linsolver(10,1); //gp,Nf
	//vd ptest={1,2,3,4};
	//a.setpResult(ptest,0);
	//tube::printVec(a.getResult());
	//Res[a.DofsPerVar+1]=2345;
	//a.setResult(Res);
	tube::printVec(a.getV());
//	printVec(Res);
return 0;
}
