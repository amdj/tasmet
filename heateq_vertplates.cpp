#include "heateq_vertplates.h"


namespace vertplates{


heateq_vertplates::heateq_vertplates(d y0,varoperations& vop):y0(y0),vop(vop),Nf(vop.Nf),Ns(vop.Ns),A(vop),uB(vop),C(vop),rho(vop),mu(vop),Tw_over_T(vop),T(vop),rho0(1),mu0(1)
{
	TRACE(9,"heateq_vertplates::heateq_vertplates(d y0,varoperations vop)");
	TRACE(9,"Ns:" << A.Ns);
	s=vdzeros(Ns+1);


}

vc heateq_vertplates::hnu_n(d y) {

	return cosh(sqrt(I)*s*y/y0)/cosh(sqrt(I)*s);
}

var heateq_vertplates::g_n(d y) {
	TRACE(9,"heateq_vertplates::g_n(d y)");
	vc hnu=hnu_n(y);
	vc gn(hnu.size());
	vc Cn=C.getcRes();
	vc Cny0=y0*Cn;
	DEBUGLOG("Cn*y0(0): "<< Cny0(0) );
	DEBUGLOG("y0:"<<y0);
	// Later on, the velocity terms have to be added...
	gn=y0*Cn%cos(y*Cn)/sin(y0*Cn);
	DEBUGLOG("y/y0:"<< y/y0);
	DEBUGLOG("gn(0):"<< gn(0));
	var gnvar(vop);
	gnvar.set(gn);
	return gnvar;
}
var heateq_vertplates::Temp(d y) {
	TRACE(9,"heateq_vertplates::Temp(d y)");
	var gn=g_n(y);
//	TRACE(9,"gn:" << gn);
//	TRACE(9,"gn timedata:" << gn.tdata());
//	TRACE(9,"T:" << T);
	var Temp=gn;
	TRACE(9,"Temp:" << Temp);
	TRACE(9,"heateq_vertplates::Temp(d y) done.");
	return Temp;
}

vd heateq_vertplates::Temp(d t,vd y) // Compute the temperature for a given time
{
	TRACE(9,"heateq_vertplates::Temp(d t,vd y)");
	us N=y.size();
	vd temp(N);
	var Ty(vop);
	d yi;
	for(us i=0; i<N;i++){
		yi=y(i);
		Ty=Temp(yi);
		temp(i)=Ty.tdata(t);
	}

	TRACE(9,"heateq_vertplates::Temp(d t,vd y) done.");
	return temp;

}

void heateq_vertplates::setdata(const var& Tw,const var& T,const var& p,const var& dTdx,const var& dpdx,gas& g)
{
	TRACE(9,"heateq_vertplates::setdata()");
	TRACE(9,"T.tdata(): " << T.tdata());
	this->T=T;
	TRACE(9,"T:" << this->T);
	Tw_over_T=Tw/T;
	var rho=g.rho(T,p);
	var mu=g.mu(T);
	var kappa=g.kappa(T);
	var gamma=g.gamma(T);


	//Compute constants to remember as class variables
	s=y0*sqrt(rho0*vop.omgvec/mu0);
	mu0=mu(0);
	rho0=rho(0);

	TRACE(9,"rho0:"<<rho0);
	TRACE(9,"mu0:" << mu0);


	vd dpdt=p.ddt().tdata();
	vd pt=p.tdata();
	vd Tt=T.tdata();
	vd dTdt=(T.ddt()).tdata();
	vd gammat=gamma.tdata();
	vd kappat=kappa.tdata();

	TRACE(9,"Kappa: " << kappa);
	vd t0=(1.0/(kappat%Tt));
	vd t1=(gammat)/(gammat-1.0)%pt/Tt;
	vd A_td=t0%(t1%dTdt-dpdt);
	vd B_td=t0%(t1%dTdx.tdata()-dpdx.tdata());
	TRACE(9,"A_td computed:" << A_td);
	TRACE(9,"B_td computed:" << B_td);

	A.settdata(A_td);
	uB.settdata(B_td);
	TRACE(9,"settdata done.");

	vd Cguess=ones<vd>(Ns)*0.55*pi/y0;
	Cguess(0)=1/y0;
	boost::function<vd (const vd&)> fun;
	fun=boost::bind( &heateq_vertplates::Cerr,this,_1);
		//fun=&(*this)::Cerr;
	vd Csol=math_anne::fsolve(Cguess,fun);
	C.set(Csol);
//	cout << "Csol:" << Csol <<endl;
	//vd Curerr=Cerr(Cguess);
	DEBUGLOG("Solution for Cn found: " << Csol);

}

vd heateq_vertplates::Cerr(const vd& C)
{
	var Cvar(vop);
	Cvar.set(C);
	vc Cn=Cvar.getcRes();
	vc An=A.getcRes();
	TRACE(9,"Cerr called with argument" << C);
	//Compute the error in C_n for current choice of C
	TRACE(9,"heateq_vertplates::Cerr(const vd& C)");
	TRACE(9,"rho0" << rho0);
	TRACE(9,"mu0" << mu0);
	TRACE(9,"Shear wave number:" << s);

	vc Tw_over_T_n=(Tw_over_T).getcRes();

	//This implementation is still without velocity
	vc Cnsq=pow(Cn,2);
	vc err_n=(Cn*y0-y0*An_over_Cnsq(An,Cn))/tan(Cn*y0)+(An_over_Cnsq(An,Cn)-Tw_over_T_n);
	var error(vop); error.set(err_n);
	return error();
}
vc heateq_vertplates::An_over_Cnsq(vc& An,vc& Cn) const
{
	vc result=vczeros(An.size());
	for (us i=0;i<Nf+1;i++){
		if(abs(An(i))<1e-12) result(i)=0;
		else if(!(Cn(i)==0.0)) {
			result(i)=An(i)/pow(Cn(i),2);
		}
		else {
			LOG4CPLUS_DEBUG(logger,"Warning: Cn(" << i << ") equals zero. This results in a singularity. Result of division artificially set to zero");
			result(i)=0;
		}
	}
	return result;

}
heateq_vertplates::~heateq_vertplates()
{
	//dtor
}
} // namespace vertplates
