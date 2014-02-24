#ifndef LINTUBE_H
#define LINTUBE_H

#include "tube.h"
#include "vvar.h"
#include "vtypes.h"
#include "freqoperators.h"
#include "continuityeq.h"
namespace tube {

class lintube : public tube {

public:
    double L,Sf;
    lintube(us gp,us Nf,double L,double Sf,double freq);
    virtual ~lintube();
	freqoperators::freqoperators fop;

	vd getx();
	void createsys();
//	void relax(us);
//	void relaxpunew(us i,vd& p,vd& u);
//		void solvesys();
protected:
	vd Result,Error;
	d K1,K2;
	dmat K;
	vd R;
	eq::eqsystem Ktot;
	eq::zeroeq zeroeqrho,zeroeqT,zeroeqTs,zeroeqp,zeroeqU;
	eq::lincontinuityeq peq;
	eq::linmomentumeq Ueq;

	us TotalDofs;
	vd x;
	double freq,omg;
	const us NumVars=2; //Number of variables involved in problem
	const d c0=343;
	const d c0sq=pow(c0,2);
	const d up=1;
	const d rho0=1.2;
	const d pm=101325;
	vd& getRes();
private:
	//eq::linmomentumeq Ueq;

//	void subsys(int i,dmat& Ksub,vd& Rsub);

//	void DLp_dpi(us i,dmat &);
//	void DLp_dui(us i,dmat &);
//	void DLp_duip1(us i,dmat&);
//	void DLp_duip2(us i,dmat&);
//	void DLp_duim2(us i,dmat&);
//	void DLp_duim1(us i,dmat&);
//
//
//	void DLu_dui(us i,dmat &);
//	void DLu_dpi(us i,dmat &);
//	void DLu_dpip1(us i,dmat&);
//	void DLu_dpip2(us i,dmat&);
//	void DLu_dpim1(us i,dmat&);
//	void DLu_dpim2(us i,dmat&);
//
//	void fillsubblock(us,us,dmat&,dmat&);
//
//	void Ksys();
//
//	void relaxp(us i,vd& dp);
//	void relaxu(us i,vd& du);
//
//	void relaxpu_firstorder(us i,vd& dp,vd& du);
};

}

#endif // LINTUBE_H
