//#define ANNE_WARNINGS 1
//#define ANNE_DEBUG_MATH 1

#include "math_anne.h"
#include "solid.h"
#include "interpolate.h"
#include "vtypes.h"
#include "reg_fast.h"

void plot(string name,vd& x,vd& y){
	std::ofstream file1;
	std::stringstream cmd1;
	cmd1 << "rm " << name << ".dat";
	system(cmd1.str().c_str());
	std::stringstream filename;
	filename << name << ".dat";
	file1.open(filename.str());
	for(us i=0;i<x.size();i++){
		file1 << x(i) << " " << y(i) << endl;
	}

	file1.close();
	std::stringstream graph;
	graph << "graph -T ps < " << filename.str() << " > " << name  << ".eps";
	//cout << graph.str();
	system(graph.str().c_str());
}


using namespace math_anne;

int main(){

//testeq testje;
//vd xguess=vdzeros(2);
//testje.bla(xguess);
////vdfunctorT<testeq> ftor(&testje,&testeq::bla);
//vdfunctorT<testeq> fun(&testje,&testeq::bla);
//eqsolver solver(&fun);
//vd x=solver(xguess);
//prn("solution:",x);
//math_anne::interpolate gc("gctable.dat");
vd x=arma::exp10(linspace<vd>(-2,2,300));
vd y(x.size());
//vd y(x.size());
for (us i=0;i<x.size();i++){ y(i)=reg::gcfun(x(i));}
plot("gc",x,y);
return 0;
}
	//vd xx=linspace<vd>(0,20.0,100);
	//vd yy=pow(xx,2);
//	math_anne::interpolate gc("gctable.dat");
	//math_anne::interpolate gv=math_anne::interpolate("gvtable.dat");


//	prn("x:",x);
//	prn("gc:",gc.call(0.5) );

//	math_anne::interpolate p1(xx,yy);
//	prn("test:",p1.call(3));
//--------------------------------------------------
//	d rh0=0.02820948;
//	d dnu0=0.000218961389962799;
//	vd rh=rh0*ones<vd>(10);
//	vc dnu=dnu0*ones<vc>(10);
//	//cout << rh/dnu << endl;
//	cout << f_blapprox(rh,dnu) << endl;
//solid::solid m=solid::solid("copper");
//vd T=linspace<vd>(200,800,10);
//cout << "Heat capacity:" << endl;
//cout << m.cs(T) << endl;
//cout << "Thermal conductivity:" << endl;
//cout << m.kappa(T) << endl;
//cout << "Density:" << endl;
//cout << m.rho(T) << endl;
//vc s=linspace<vc>(-2,-2,1).transform([](c x){return pow(10,x);});
//
//pyprn("s:",s);
//vd Pr=0.68*ones<vd>(1);
//pyprn("Pr:",Pr);
//vc eps_s=zeros<vc>(1);
//math_anne::ThetaLambda t=math_anne::ThetaLambda();
//vc thet=t.Theta(s,Pr,eps_s);
//vd iml=t.ImLambda(s,Pr,eps_s);
//pyprn("Theta:",thet);
//pyprn("ImLambda:",iml);
//pyprn("s=0.1,Pr=0.68 imag Theta:",imag(t.Theta(s,Pr,eps_s)));
