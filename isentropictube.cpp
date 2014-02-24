#include "isentropictube.h"

namespace tube{

isentropictube::isentropictube(us gp1,us Nf1,double L1,double Sf1,double freq1,double Up1): tube(gp1,Nf1),L(L1),Sf(Sf1),freq(freq1),Up(Up1) {
    us i;
	omg=2*pi*freq1;
	double h=L/double(gp);
	xp=var::var(Nf);
	pp=var::var(Nf,0);
	xp.setRes(-Up/omg,2);
	//assert((Up/omg)<(h/2));
	V0=var::var(Nf,h*Sf);
	V0.setRes(V0.getTdata()-xp.getTdata()*Sf);

    for(i=0; i<gp; i++) {
		// Not exactly on the gridpoints, but inbetween
		x(i)=0.5*h+double(i)*h;
		V(i)=h;
//		Utest[i]=Up*cos(2*pi*x[i]/L);
    }


	TotalDofs=2*DofsPerVar; //For the lintube case, since we only require pressure and velocity.
	//cout << "TotalDofs:" << TotalDofs << endl;
	Result=vd(TotalDofs);
	Error=vd(TotalDofs);
}
vd isentropictube::getpError(us i) {
	vd pError=vd(Ns);
	if(i>0 && i<gp-1) { //Interior node
		vd drhodt=rho[i].ddt_fd(omg);
		pError=V[i]*drhodt+rhoU.dotn(i,x);

	}
	else if(i==0) {//Treat first node
		var::var rhoV=var::var(Nf);
		rhoV.setTdata(rho[0].getTdata()*V0.getTdata());
		vd rhoUdotn=0.5*(rhoU[i].getAdata()+rhoU[i+1].getAdata());
		pError=rhoV.ddt_fd(omg)+rhoUdotn;
		pError[0]=p[0].getAdata(0)-pm;
	}
	else {//Right node
		vd drhodt=rho[i].ddt_fd(omg);
		vd rhoUdotn=-0.5*(rhoU[i].getAdata()+rhoU[i-1].getAdata());
		pError=V[i]*drhodt+rhoUdotn;
	}

	////////////LINTUBEFV SHIT START
//		vd dpdt=p[i].ddt_fd(omg);
//	vd Ui=U[i].getAdata();
//	if(i>0 &&i<gp-1) {
//		vd Udotn=U.dotn(i,x);
//		pError=V[i]*dpdt/c0sq+rho0*Udotn;
//	}
//	else if(i==0) { //First node
//		vd Uip1=U[i+1].getAdata();
//		vd Udotn=0.5*(Uip1+Ui);
//		pError=V[i]*dpdt/c0sq+rho0*V0.ddt_fd(omg)+rho0*Udotn;
//	}
//	else //Right node
//	{
//		vd Uim1=U[i-1].getAdata();
//		vd Udotn=-0.5*(Uim1+Ui);
//		pError=V[i]*dpdt/c0sq+rho0*Udotn;
//	}
//	//if(i==gp-1){ cout << "i=gp-1" << endl;
	pError[0]=p[i].getAdata()[0]-101325;
	////////////LINTUBEFV SHIT END
	return pError;
}
vd isentropictube::getUError(us i) {
	vd UError=vd(Ns);
	vd rhoUudotn(Ns);
	if(i>0 && i<gp-1) { //Interior node
		vd drhoUdt=rhoU[i].ddt_fd(omg);
		UError=V[i]*drhoUdt+rhoUu.dotn(i,x)+p.dotn(i,x); //Should be the right one
		//UError=V[i]*drhoUdt+p.dotn(i,x); // For testing, we leave convective term out
	}
	else if(i==0) {//Treat first node, special piston problem
		var::var rhoUV=var::var(Nf);
		rhoUV.setTdata((rhoU[0].getTdata())*(V0.getTdata()));
		double h=x[1]-x[0];
		vd phalf=0.5*(p[0].getAdata()+p[1].getAdata());
		vd phalf_td=0.5*(p[0].getTdata()+p[1].getTdata());
		vd dpdxhalf_td=(p[1].getTdata()-p[0].getTdata())/h;
		vd hvec(Ns);
		us k;
		for(k=0;k<Ns;k++) { hvec[k]=h;}
		vd h_till_piston=hvec-xp.getTdata();
		pp.setTdata(phalf_td-dpdxhalf_td*h_till_piston);
		vd pL=pp.getAdata();
		rhoUudotn=0.5*(rhoUu[0].getAdata()+rhoUu[1].getAdata());
		vd pdotn=phalf-pL;
		UError=rhoUV.ddt_fd(omg)+rhoUudotn+pdotn;

	}
	else {//Right node
		vd drhoUdt=rhoU[i].ddt_fd(omg);
		rhoUudotn=-0.5*(rhoUu[i-1].getAdata()+rhoUu[i].getAdata());
		UError=V[i]*drhoUdt+rhoUudotn+p.dotn(i,x);
//		UError=drhoUdt+p.dotn(i,x);// For testing, we leave convective term out
	}
	UError[0]=U[i].getAdata(0); //Explicitly put time-avg velocity to zero, if not Jacobian becomes zero!
	//Fill in shit of lintubelv
	//vd pdotn=p.dotn(i,x);
	//vd dUdt=U[i].ddt_fd(omg);
	//UError=rho0*dUdt*V[i]+pdotn;
	UError[0]=U[i].getAdata(0); //Explicitly put time-avg velocity to zero, if not Jacobian becomes zero!
	//End of shit of lintubelv
	return UError;
}
vd isentropictube::getUError() {
	us i;
	for(i=0;i<gp;i++)	U.Error.subvec(i*Ns,(i+1)*Ns)=getUError(i);
	return U.Error;
}
vd isentropictube::getpError() {
	us i;
	for(i=0;i<gp;i++)	p.Error.subvec(i*Ns,(i+1)*Ns)=getpError(i);
	return p.Error;
}
void isentropictube::setResult(vd Vec) {
	U.setResult(Vec.subvec(0,DofsPerVar));
	p.setResult(Vec.subvec(DofsPerVar,2*DofsPerVar));
	rho.setTdata(rho0*pow(p.getTdata()/pm,1/gamm)); //Something wrong here???
	rhoU.setTdata(rho.getTdata()%U.getTdata());
	rhoUu.setTdata(rhoU.getTdata()%U.getTdata()/Sf);
}
vd isentropictube::getResult() {
	Result.subvec(0,DofsPerVar)=U.getResult();
	Result.subvec(DofsPerVar,2*DofsPerVar)=p.getResult();
	return Result;
}
vd isentropictube::getError() {
	Error.subvec(0,DofsPerVar)=getUError();
	Error.subvec(DofsPerVar,2*DofsPerVar)=getpError();
	return Error;
}
vd isentropictube::getx()  {return x; }
isentropictube::~isentropictube() {}

}
