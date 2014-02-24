#include "continuityeq.h"

namespace eq{

subsys::subsys(){}
subsys::subsys(us Nf): Nf(Nf) {
	Ns=2*Nf+1;
	bw=Ns-1;
	reset();
	//Zeroth: density
	//First : Volume flow
	//Second: Temperature
	//Third:  Pressure
	//Fourth: Solid temperature

}
void subsys::setR(vd& data){
	Rsub=data;
}
void subsys::reset(){
	Ksub=zeros<dmat>(Ns,3*NVARS*Ns);
	Rsub=zeros<vd>(Ns);
}
dmat subsys::subblock(us pos){
	return Ksub.submat(0,pos*Ns,bw,pos*Ns+bw);
}
void subsys::setsubblock(us pos,vd& data){
	Ksub.submat(0,pos*Ns,bw,pos*Ns+bw).diag()=data;
}
void subsys::setsubblock(us pos,dmat& data){
	Ksub.submat(0,pos*Ns,bw,pos*Ns+bw)=data;
}
dmat subsys::rhoblock(us i){
	return subblock((i+1)*NVARS+0);
}
dmat subsys::Ublock(us i){
	return subblock((i+1)*NVARS+1);
}
dmat subsys::Tblock(us i){
	return subblock((i+1)*NVARS+2);
}
dmat subsys::pblock(us i){
	return subblock((i+1)*NVARS+3);
}
dmat subsys::Tsblock(us i){
	return subblock((i+1)*NVARS+4);
}
void subsys::setrhoblock(us i,vd& data){
	setsubblock(NVARS*(i+1),data);
}
void subsys::setUblock(us i,vd& data){
	setsubblock(NVARS*(i+1)+1,data);
}
void subsys::setTblock(us i,vd& data){
	setsubblock(NVARS*(i+1)+2,data);
}
void subsys::setpblock(us i,vd& data){
	setsubblock(NVARS*(i+1)+3,data);
}
void subsys::setTsblock(us i,vd& data){
	setsubblock(NVARS*(i+1)+4,data);
}
void subsys::setrhoblock(us i,dmat& data){
	setsubblock(NVARS*(i+1),data);
}
void subsys::setUblock(us i,dmat& data){
	setsubblock(NVARS*(i+1)+1,data);
}
void subsys::setTblock(us i,dmat& data){
	setsubblock(NVARS*(i+1)+2,data);
}
void subsys::setpblock(us i,dmat& data){
	setsubblock(NVARS*(i+1)+3,data);
}
void subsys::setTsblock(us i,dmat& data){
	setsubblock(NVARS*(i+1)+4,data);
}
dmat subsys::getK() {return Ksub;}
vd subsys::getR() {return Rsub;}
subsys::~subsys(){}

eqsystem::eqsystem(){}
eqsystem::eqsystem(us Nf,us gp): Nf(Nf),gp(gp){
	Ns=2*Nf+1;
	bh=Ns-1; //Block height
	bw=3*NVARS*Ns-1; //Block width
	reset();
}
void eqsystem::reset(){

	K=zeros<dmat>(NVARS*Ns*gp,NVARS*Ns*gp);
	R=zeros<vd>(NVARS*Ns*gp);
}
dmat& eqsystem::getK(){return K;}
vd& eqsystem::getR(){return R;}
void eqsystem::seteq(us eqnr,us gpnr,subsys& Ksub){
	us firstrow;
	us firstcol;
	if(gpnr==0){//First node
		firstrow=(eqnr)*Ns;
		firstcol=0;
	}
	else if(gpnr>0&&gpnr<gp-1){
//		prn("Middle node",0);
//		prn("gpnr:",gpnr);
//		prn("eqnr:",eqnr);
		firstrow=(gpnr*NVARS+eqnr)*Ns;
		firstcol=((gpnr-1)*NVARS)*Ns;
//		prn("Firtrow:",firstrow);
//		prn("Firtcol:",firstcol);
	}
	else{ //Last node

		firstrow=(gpnr*NVARS+eqnr)*Ns;
		firstcol=((gpnr-2)*NVARS)*Ns;
	}
	us lastrow=firstrow+bh;
	us lastcol=firstcol+bw;
	K.submat(firstrow,firstcol,lastrow,lastcol)=Ksub.getK();
//	prn("Ksub:",Ksub.getK());
	R.subvec(firstrow,lastrow)=Ksub.getR();

}
eqsystem::~eqsystem(){}

equation::equation(){}
equation::equation(us gp,us Nf,d freq):
 gp(gp),freqoperators(Nf,freq) {
	bw=Ns-1; //Blockwidth
	Ksub=subsys(Nf);
}

zeroeq::zeroeq(){}
zeroeq::zeroeq(us gp,us Nf,us varnr): equation(gp,Nf,1.0),varnr(varnr){}
subsys zeroeq::getsys(us gpnr){
	Ksub.reset();
	dmat identity=eye<dmat>(Ns,Ns);
	if(gpnr==0){
		Ksub.setsubblock(varnr,identity);
	}
	else if(gpnr>0 &&gpnr<gp-1){
		Ksub.setsubblock(NVARS+varnr,identity);
	}
	else {
		Ksub.setsubblock(2*NVARS+varnr,identity);
	}
	return Ksub;
}

lincontinuityeq::lincontinuityeq(){

}
lincontinuityeq::lincontinuityeq(us gp,us Nf,vd& x,d freq):equation(gp,Nf,freq),x(x) {

    d dx=x[2]-x[1];
    K1=rho0*c0sq/2/dx;
    K2=1/(rho0*dx);
}

subsys lincontinuityeq::getsys(us gpnr){
	//Create the local matrix
	Ksub.reset();
	dmat K1m=-K1*eye<dmat>(Ns,Ns);
	dmat K1p=K1*eye<dmat>(Ns,Ns);
	vd ZEROS=vdzeros(Ns);
	vd Rlid=ZEROS;
	dmat zerosrow=zeros<dmat>(1,Ns);
	if(gpnr>0 &&gpnr<gp-1){

		dmat pblk=DDTfd;
		pblk.row(0)=zerosrow;
		pblk(0,0)=1.0;

		Rlid(0)=pm;
		Ksub.setUblock(-1,K1m);
		Ksub.setpblock(0,pblk);
		Ksub.setUblock(1,K1p);
		Ksub.setR(Rlid);
	}
	else if(gpnr==0){
		dmat K10=-3*K1p;
		dmat K11=4*K1p;
		dmat K12=-K1p;

		vd ZEROS=vdzeros(Ns);

		K10.row(0)=zerosrow;
		K11.row(0)=zerosrow;
		K12.row(0)=zerosrow;
		Rlid(0)=pm;

		dmat pblk=DDTfd;
		pblk.row(0)=zerosrow;
		pblk(0,0)=1.0;

		Ksub.setpblock(-1,pblk);
		Ksub.setUblock(-1,K10);
		Ksub.setUblock(0,K11);
		Ksub.setUblock(1,K12);
		Ksub.setR(Rlid);
	}
	else {
		dmat KN0=3*K1p;
		dmat KNm1=-4*K1p;
		dmat KNm2=K1p;
		Ksub.setUblock(-1,KNm2);
		Ksub.setUblock(0,KNm1);
		Ksub.setpblock(1,DDTfd);
		Ksub.setUblock(1,KN0);
	}
	return Ksub;
}
lincontinuityeq::~lincontinuityeq(){
	//dtor
}

linmomentumeq::linmomentumeq(){

}
linmomentumeq::linmomentumeq(us gp,us Nf,vd& x,d freq):equation(gp,Nf,freq),x(x) {

    d dx=x[2]-x[1];
    K1=rho0*c0sq/2/dx;
    K2=1/(2*rho0*dx);
}
subsys linmomentumeq::getsys(us gpnr){
	//Create the local matrix
	Ksub.reset();
	dmat K2m=-K2*eye<dmat>(Ns,Ns);
	dmat K2p=K2*eye<dmat>(Ns,Ns);
	dmat EYE=eye<dmat>(Ns,Ns);
	if(gpnr>0 &&gpnr<gp-1){
		Ksub.setpblock(-1,K2m);
		Ksub.setUblock(0,DDTfd);
		Ksub.setpblock(1,K2p);
	}
	else if(gpnr==0){

		Ksub.setUblock(-1,EYE);
		vd leftbc=vdzeros(Ns); leftbc(1)=1;
		Ksub.setR(leftbc);
	}
	else {
		Ksub.setUblock(1,EYE);
		vd rightbc=vdzeros(Ns); rightbc(1)=0;
		Ksub.setR(rightbc);
	}
	return Ksub;
}
linmomentumeq::~linmomentumeq(){
	//dtor
}

solution::solution(){}
solution::solution(us gp,us Nf,vd& sol):Nf(Nf),gp(gp){
	Ns=2*Nf+1;
	Usol=vdzeros(Ns*gp);
	psol=vdzeros(Ns*gp);
	setsol(sol);
}
void solution::setsol(vd& sol1){
	sol=sol1;
	Usol=getvarsol(Unr);
	psol=getvarsol(pnr);
}
vd solution::getvarsol(us varnr){
	us firstrow;
	us lastrow;
	vd varsol=vdzeros(gp*Ns);
	for(us i=0;i<gp;i++){
		firstrow=(i*NVARS+varnr)*Ns;
		lastrow=firstrow+Ns-1;
		varsol.subvec(i*Ns,(i+1)*Ns-1)=sol.subvec(firstrow,lastrow);
	}
	return varsol;
}
vd solution::getU(us freqnr){
	vd Usolfreq=vdzeros(gp);
	for(us i=0;i<gp;i++){
		Usolfreq(i)=Usol(Ns*i+freqnr);
	}
	return Usolfreq;
}
vd solution::getp(us freqnr){
	vd psolfreq=vdzeros(gp);
	for(us i=0;i<gp;i++){
		psolfreq(i)=psol(Ns*i+freqnr);
	}
	return psolfreq;
}
} //namespace eq
