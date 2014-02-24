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


//void lintube::relaxpu_firstorder(us i,vd& dp,vd& du){
//	dmat Ksys=zeros<dmat>(2*(Ns-1),2*(Ns-1));
//	vd R1=vd(Ns-1);
//	vd R2=vd(Ns-1);
//	vd pold=p[i].getResfluc();
//	vd uold=u[i].getResfluc();
//
//	if (i==0){
//		d dx=x[1]-x[0];
//		d K1=rho0*c0sq/(dx);
//		d K2=1.0/(rho0*dx);
//		vd u0=vd(Ns-1); u0(0)=up;
//		du=u0-uold;
//		vd uright=u[i+1].getResfluc();
//		dp=K1*(uright-uold)-fop.ddt*pold;
//	}
//	else if (i>0&&i<gp-1){
//		d dx=x[i]-x[i-1];
//		d K1=rho0*c0sq/(dx);
//		d K2=1.0/(rho0*dx);
//
//
//		vd uleft=u[i-1].getResfluc();
//		vd uright=u[i+1].getResfluc();
//		vd pleft=p[i-1].getResfluc();
//		vd pright=p[i+1].getResfluc();
//
//
//		R1=K1*uleft-fop.ddt*pold;
//		R2=K2*pleft-fop.ddt*uold;
//		vd R=vd(2*(Ns-1));
//
//
//		Ksys.submat(0,0,Ns-2,Ns-2)=fop.ddt;
//		Ksys.submat(Ns-1,Ns-1,2*(Ns-1)-1,2*(Ns-1)-1)=fop.ddt;
//		Ksys.submat(0,Ns-1 , Ns-2, 2*(Ns-1)-1).diag()+=K1;
//		Ksys.submat(Ns-1,0, 2*(Ns-1)-1,Ns-2).diag()+=K2;
//		R.subvec(0,Ns-2)=R1;
//		R.subvec(Ns-1,2*(Ns-1)-1)=R2;
//		//cout << Ksys;
//
//		vd dpdu=solve(Ksys,R);
//		//cout << dpdu << endl;
//		dp=dpdu.subvec(0,Ns-2);
//		du=dpdu.subvec(Ns-1,2*(Ns-1)-1);
//	}
//	else if(i==gp-1){
//		d dx=x[i]-x[i-1];
//		d K1=rho0*c0sq/(dx);
//		d K2=1.0/(rho0*dx);
//		vd u0=vd(Ns-1);
//		du=u0-uold;
//		vd uleft=u[i-1].getResfluc();
//		dp=K1*(uold-uleft)-fop.ddt*pold;
//	}
//
//}
//void lintube::relax(us number){
//	vd dp=vd(Ns-1);
//	vd du=vd(Ns-1);
//	vd pold=vd(Ns-1);
//	vd uold=vd(Ns-1);
//	vd pnew=vd(Ns-1);
//	vd unew=vd(Ns-1);
//	vd pleft=vd(Ns-1);
//	vd pright=vd(Ns-1);
//	vd uleft=vd(Ns-1);
//	vd uright=vd(Ns-1);
//	vd pleftnew=vd(Ns-1);
//	vd prightnew=vd(Ns-1);
//	vd uleftnew=vd(Ns-1);
//	vd urightnew=vd(Ns-1);
//
//	for (us n=0;n<number;n++){
//
//		//us i=0;
//		for (us i=0;i<gp;i++){
//			pold=p[i].getResfluc();
//			uold=u[i].getResfluc();
//			relaxpu_firstorder(i,dp,du);
//			pnew=pold+dp;
//			unew=uold+du;
//			p[i].setResfluc(pnew);
//			//u[i].setResfluc(unew);
//		}
////
////		relaxp(i,dp);
////		relaxu(i,du);
////		pnew=pold+dp;
////		unew=uold+du;
////
//////		p[i].setResfluc(pnew);
////		u[i].setResfluc(unew);
////
////		for(i=1;i<gp-1;i++){
////			uleft=u[i-1].getResfluc();
////			uright=u[i+1].getResfluc();
////			pleft=p[i-1].getResfluc();
////			pright=p[i+1].getResfluc();
////			pold=p[i].getResfluc();
////			uold=u[i].getResfluc();
//////			relaxpunew(i,dp,du);
//////
////
////			relaxp(i,dp);
////			relaxu(i,du);
////			pnew=pold+0.6*dp;
////			unew=uold+0.6*du;
////
////			pleftnew=pleft-dp/2;
////			prightnew=pright+dp/2;
////			uleftnew=uleft-du/2;
////			urightnew=uright+du/2;
////
////
////			u[i].setResfluc(unew);
////			p[i].setResfluc(pnew);
////			u[i-1].setResfluc(uleftnew);
////			p[i-1].setResfluc(pleftnew);
////			u[i+1].setResfluc(urightnew);
////			p[i+1].setResfluc(prightnew);
////		}
//////
////		i=gp-1;
////		pold=p[i].getResfluc();
////		uold=u[i].getResfluc();
////		relaxp(i,dp);
////		relaxu(i,du);
////		pnew=pold+dp;
////		unew=uold+du;
//////		p[i].setResfluc(pnew);
////		u[i].setResfluc(unew);
//
////
//
////		}
//	}
//
//
//}
//void lintube::relaxp(us i,vd& dp){ //Update density in node i-1 and node i such that continuity equation is satisfied
//
////	vd dp=vd(Ns-1); //One smaller, because steady component excluded
//	p[i].setRes(pm,0); //Update steady component.
//
////	if (i>0 && i<gp-1){ //interior nodes
//
////		vd poldip1=p[i+1].getRes().subvec(1,Ns-1); //Only extract fluctuating components
////		vd poldim1=p[i-1].getRes().subvec(1,Ns-1); //Only extract fluctuating components
////		vd uiold=u[i].getRes().subvec(1,Ns-1);
//		//cout << uiold;
////		d twodx=x[i+1]-x[i-1];
//
////		vd dp = 0.5*(twodx*rho0*fop.ddt*uiold-poldip1+poldim1);
////		vd pnewip1=poldip1+dp;
//////		vd pnewim1=poldip1-dp;
////		p[i+1].setResfluc(pnewip1);
////		p[i-1].setResfluc(pnewim1);
//
////	}
////	else if	(i==0) { //Left boundary
//
////	}
////	else if(i==1) {}
////	cout << "pold:" << pold << endl;
////	cout << "ddt:" << fop.ddt << endl;
//	//cout << "pold:" << pold << endl;
//	vd pold=p[i].getResfluc();
////	vd dp=vd(Ns-1);
//	vd dudx=u.ddx_central(i,x).subvec(1,Ns-1);
//	dp = fop.iddt*(-1.0*rho0*c0sq*dudx-fop.ddt*pold);
////	vd pnew=pold+relaxomega*dp;
//	//cout << "pnew:" << pnew << endl;
////	p[i].setResfluc(pnew);
//}
//void lintube::relaxu(us i,vd& du){
//	u[i].setRes(0,0.0); //Set result for steady component to zero;
//	vd uold=u[i].getResfluc(); //Only extract
//	if(i>0 && i<gp-1){
//		//vd du=vd(Ns-1);
//		//cout << "uold:" << uold << endl;
//		vd dpdx=p.ddx_backward(i,x).subvec(1,Ns-1);
//		//vd uleft=u[i-1].getRes().subvec(1,Ns-1);
//		du=fop.iddt*(-1.0*dpdx/rho0-fop.ddt*uold);
//		//du=uleft-uold;
////		vd unew=uold+relaxomega*du;
//		//cout << "unew:" << unew << endl;
////		u[i].setResfluc(unew);
//	}
//	else if(i==0){
//		vd u0=vd(Ns-1); u0(0)=up;
//		du=u0-uold;
////		u[i].setResfluc(u0);
//	}
//	else if(i==gp-1){
//		vd u0=vd(Ns-1);
//		du=u0-uold;
////		u[i].setResfluc(u0);
//	}
//}
//
//
//void lintube::solvesys(){
//	dmat Ksys=zeros<dmat>(TotalDofs);
//	dmat Ksub=zeros<dmat>(2*Ns,6*Ns);
//	vd Rsub=zeros<vd>(2*Ns);
//
//	subsys(1,Ksub,Rsub);
//	cout << Ksub << endl;
//
//}

//void lintube::Ksys(){
//	K=zeros<dmat>(TotalDofs,TotalDofs);
//	dmat Ksub=zeros<dmat>(Ns,Ns);
//	vd Rsub=zeros<vd>(Ns);
//	for(us i=0;i<gp;i++){
//		subsysInner(i, Ksub, Rsub);
//
//	}
//
//}


//void lintube::subsys(int i,dmat& Ksub,vd& Rsub){
//	Ksub=zeros<dmat>(2*Ns,6*Ns);
//	Rsub=zeros<vd>(6*Ns);
//
//	dmat dLp_dpi=zeros<dmat>(Ns,Ns); DLp_dpi(i,dLp_dpi);
//	dmat dLp_dui=zeros<dmat>(Ns,Ns); DLp_dui(i,dLp_dui);
//	dmat dLp_duip1=zeros<dmat>(Ns,Ns); DLp_duip1(i,dLp_duip1);
//	dmat dLp_duim1=zeros<dmat>(Ns,Ns); DLp_duim1(i,dLp_duim1);
//	dmat dLp_duim2=zeros<dmat>(Ns,Ns); DLp_duim2(i,dLp_duim2);
//	dmat dLp_duip2=zeros<dmat>(Ns,Ns); DLp_duip2(i,dLp_duip2);
//
//	dmat dLu_dui=zeros<dmat>(Ns,Ns); DLu_dui(i,dLu_dui);
//	dmat dLu_dpi=zeros<dmat>(Ns,Ns); DLu_dpi(i,dLu_dpi);
//	dmat dLu_dpip1=zeros<dmat>(Ns,Ns); DLu_dpip1(i,dLu_dpip1);
//	dmat dLu_dpim1=zeros<dmat>(Ns,Ns); DLu_dpim1(i,dLu_dpim1);
//	dmat dLu_dpim2=zeros<dmat>(Ns,Ns); DLu_dpim2(i,dLu_dpim2);
//	dmat dLu_dpip2=zeros<dmat>(Ns,Ns); DLu_dpip2(i,dLu_dpip2);
//
//	if (i>0&&i<gp-1){
//		//First row, second block: dLpduim1
//		fillsubblock(1,2,Ksub,dLp_duim1);
//		//First row, third block: dLpdpi
//		fillsubblock(1,3,Ksub,dLp_dpi);
//		//First row, sixth block: dLpduip1
//		fillsubblock(1,6,Ksub,dLp_duip1);
//
//		//Second row, first block: dLudpim1
//		fillsubblock(2,1,Ksub,dLu_dpim1);
//		//Second row, fourth block: dLu_du1
//		fillsubblock(2,4,Ksub,dLu_dui);
//		//Second row, fifth block: dLu_dpip1
//		fillsubblock(2,5,Ksub,dLu_dpip1);
//	}
//	else if(i==0){
//		//For the first node, the submatrix will be shifted.
//
//		dmat Identity=ones<dmat>(Ns,Ns);
//		fillsubblock(2,2,Ksub,Identity);
//		Rsub(Ns+1)=up;
//		fillsubblock(1,1,Ksub,dLp_dpi);
//		fillsubblock(1,2,Ksub,dLp_dui);
//		fillsubblock(1,4,Ksub,dLp_duip1);
//		fillsubblock(1,6,Ksub,dLp_duip2);
//	}
//	else if(i==gp-1){
//		dmat Identity=ones<dmat>(Ns,Ns);
//		fillsubblock(2,6,Ksub,Identity); //Set velocity to zero
//		fillsubblock(1,6,Ksub,dLp_dui);
//		fillsubblock(1,4,Ksub,dLp_duim1);
//		fillsubblock(1,2,Ksub,dLp_duim2);
//	}
//}




//void lintube::fillsubblock(us row,us col,dmat& bigmat,dmat& data) {
//	bigmat.submat((row-1)*Ns,(col-1)*Ns,(row)*Ns-1,(col)*Ns-1)=data;
//	}

//Pressure equations
//void lintube::DLp_dpi(us i,dmat& dLp_dpi){
//	dLp_dpi=fop.DDTfd;
//}
//void lintube::DLp_dui(us i,dmat& dLp_dui){
//	if(i==0) { dLp_dui=-3.0*K1*eye<dmat>(Ns,Ns); }
//	else if(i==gp-1){dLp_dui=3.0*K1*eye<dmat>(Ns,Ns); }
//	else {dLp_dui=zeros<dmat>(Ns,Ns);}
//}
//void lintube::DLp_duip1(us i,dmat& dLp_duip1){
//	assert(i<gp);
//	if (i>0&&i<gp-1){
//		dLp_duip1=K1*eye<dmat>(Ns,Ns);
//		return;
//	}
//	else if(i==0){
//		dLp_duip1=4*K1*eye<dmat>(Ns,Ns);
//		return;
//	}
//
//}
//void lintube::DLp_duim1(us i,dmat& dLp_duim1){
//	assert(i>0);
//	if (i>0&&i<gp-1){
//		dLp_duim1=-K1*eye<dmat>(Ns,Ns);
//	}
//	else if (i==gp-1){
//		dLp_duim1=3.0*K1*eye<dmat>(Ns,Ns);
//	}
//}
//void lintube::DLp_duim2(us i,dmat& dLp_duim2){
//	assert(i>0);
//	if (i==gp-1){
//		dLp_duim2=-1.0*K1*eye<dmat>(Ns,Ns);
//	}
//	else{ dLp_duim2=zeros<dmat>(Ns,Ns);}
//}
//void lintube::DLp_duip2(us i,dmat& dLp_duip2){
//	assert(i>0);
//	if (i==0){
//		dLp_duip2=-1.0*K1*eye<dmat>(Ns,Ns);
//	}
//	else{ dLp_duip2=zeros<dmat>(Ns,Ns);}
//}

//Velocity equations
//void lintube::DLu_dui(us i,dmat& dLu_dui){
//	dLu_dui=fop.DDTfd;
//}
//void lintube::DLu_dpi(us i,dmat& dLu_dpi){
//	if(i==0) { dLu_dpi=-3.0*K2*eye<dmat>(Ns,Ns); }
//	else if(i==gp-1){dLu_dpi=3.0*K2*eye<dmat>(Ns,Ns); }
//	else {dLu_dpi=zeros<dmat>(Ns,Ns);}
//}
//void lintube::DLu_dpim1(us i,dmat& dLu_dpim1){
//	assert(i>0);
//
//	if(i==gp-1){
//		dLu_dpim1=-4.0*K2*eye<dmat>(Ns,Ns);
//	}
//	else{
//		dLu_dpim1=-K2*eye<dmat>(Ns,Ns);
//	}
//}
//void lintube::DLu_dpip1(us i,dmat& dLu_dpip1){
//	assert(i<gp-1);
//	if(i==0){
//		dLu_dpip1=4.0*K2*eye<dmat>(Ns,Ns);
//	}
//	else{
//		dLu_dpip1=K2*eye<dmat>(Ns,Ns);
//	}
//}
//void lintube::DLu_dpip2(us i,dmat& dLu_dpip2){
//	assert(i>0);
//	if (i==0){
//		dLu_dpip2=-1.0*K2*eye<dmat>(Ns,Ns);
//	}
//	else{ dLu_dpip2=zeros<dmat>(Ns,Ns);}
//}
//void lintube::DLu_dpim2(us i,dmat& dLu_dpim2){
//	assert(i>0);
//	if (i==gp-1){
//		dLu_dpim2=-1.0*K2*eye<dmat>(Ns,Ns);
//	}
//	else{ dLu_dpim2=zeros<dmat>(Ns,Ns);}
//}

//vd& lintube::getRes() {
//	Result.subvec(0,DofsPerVar-1)=U.getRes();
//	Result.subvec(DofsPerVar,2*DofsPerVar)=p.getRes();
//	return Result;
//}


vd lintube::getx()  {return x; }
//vd lintube::getV()  {return V; }

lintube::~lintube() {
    //dtor
}


} // Namespace tube




