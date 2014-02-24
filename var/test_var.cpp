#include "var.h"
#include "inlinefcn.h"
#include "vvar.h"
using namespace std;


int main() {
	using var::var;
    us Nf=2;
    us gp=10;
    double f=10;
//    cout << 'hoi' << endl;
    double omg=2*pi*f;
    double T=1/f;
    //var v=var(Nf);
	vvar::vvar p=vvar::vvar("Pressure",gp,Nf,101325);
	cout << p.name() << endl;
	//printVec(p[0].getAdata());
	printVec(p.getRes());
    //double tdata[] = {11,10.6235,9.77748,9.09903,9.09903,9.77748,10.6235};

    //complex<double> I=(0,1);



    return 0;
}
