#include "gas.h"






int main(){

	initlog(log4cplus::ALL_LOG_LEVEL);


	gases::gas a("air");
	LOG4CPLUS_DEBUG(logger,"Test");
	vd T(10);
	d p=101325;
	T.fill(293.15);

	cout << a.rho(T,p) << endl;
	cout << a.cp(T) << endl;
	cout << a.mu(T) << endl;
	cout << a.kappa(T) << endl;


	cin.get();
	return 0;
}
