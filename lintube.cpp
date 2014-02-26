#include "lintube.h"
#include "continuityeq.h"
namespace tube {

lintube::lintube(us gp,us Nf,double L,double Sf,double freq)
 : tube(gp,Nf),L(L),Sf(Sf),freq(freq) {
    unsigned i;
	omg=2*pi*freq;
    double h=L/double(gp-1);
    x=vd(gp);
    for(i=0; i<gp; i++) {
		x(i)=double(i)*h;
    }

    d dx=x[2]-x[1];
    K1=rho0*c0sq/2/dx;
    K2=1/(rho0*dx);

	TotalDofs=NVARS*DofsPerVar; //For the lintube case, since we only require pressure and velocity.
	Result=vd(TotalDofs);
	Error=vd(TotalDofs);
	fop=freqoperators::freqoperators(Nf,freq);

	peq=eq::lincontinuityeq(gp,Nf,x,freq);
	Ueq=eq::linmomentumeq(gp,Nf,x,freq);
	zeroeqrho=eq::zeroeq(gp,Nf,0);

	zeroeqT=eq::zeroeq(gp,Nf,2);
	zeroeqp=eq::zeroeq(gp,Nf,3);
	zeroeqTs=eq::zeroeq(gp,Nf,4);

	Ktot=eq::eqsystem(Nf,gp);

//	prn("Ktot:",Ktot.getK().submat(0,0,NVARS*Ns,NVARS*Ns));
	//prn("Ktotrows:",Ktot.getK().n_rows);
}
void lintube::createsys(){
	eq::subsys rhosub,psub,Tsub,Tssub,Usub;
	rhosub=zeroeqrho.getsys(0);
	//prn("rhosub0:",rhosub.subblock(0));

	for(us gpnr=0;gpnr<gp;gpnr++){
		rhosub=zeroeqrho.getsys(gpnr);
		Usub=Ueq.getsys(gpnr);
		Tsub=zeroeqT.getsys(gpnr);
		psub=peq.getsys(gpnr);
		Tssub=zeroeqTs.getsys(gpnr);
		Ktot.seteq(0,gpnr,rhosub);
		Ktot.seteq(1,gpnr,Usub);
		Ktot.seteq(2,gpnr,Tsub);
		Ktot.seteq(3,gpnr,psub);
		Ktot.seteq(4,gpnr,Tssub);
	}
	dmat K=Ktot.getK();

//	prn("Diag k:",Ktot.getK().diag());
//	prn("Diag-1 K:",Ktot.getK().diag(-1));
//	prn("Diag+1 K:",Ktot.getK().diag(1));


	//prn("Det k:",detk);
	//K.diag()=zeros<vd>(TotalDofs);
	cout << "Starting to solve the system..." << endl;
	vd sol=solve(Ktot.getK(),Ktot.getR());
	//prn("Sol:",sol);



	eq::solution solu=eq::solution(gp,Nf,sol);

	pyprn("U0:",solu.getU(0));
	pyprn("U1cos:",solu.getU(1));
	//pyprn("U1sin:",solu.getU(2));
	//pyprn("p0:",solu.getp(0));
	//pyprn("p1cos:",solu.getp(1));
	pyprn("p1sin:",solu.getp(2));
	//prn("K:",K);
	//prn("Ktot:",norm(K,2));
}



vd lintube::getx()  {return x; }
//vd lintube::getV()  {return V; }

lintube::~lintube() {
    //dtor
}


} // Namespace tube




