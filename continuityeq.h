#ifndef CONTINUITYEQ_H
#define CONTINUITYEQ_H
#include "vtypes.h"
#include "freqoperators.h"

#define NVARS 5

namespace eq{

class subsys{
	public:
		subsys();
		subsys(us Nf);
		~subsys();
		void setR(vd& data);
		dmat rhoblock(us);
		dmat Ublock(us);
		dmat Tblock(us);
		dmat pblock(us);
		dmat Tsblock(us);
		void setrhoblock(us,vd&);
		void setUblock(us,vd& );
		void setTblock(us,vd& );
		void setpblock(us,vd& );
		void setTsblock(us,vd& );
		void setrhoblock(us,dmat&);
		void setUblock(us,dmat& );
		void setTblock(us,dmat& );
		void setpblock(us,dmat& );
		void setTsblock(us,dmat& );
		void reset();
		dmat getK();
		vd getR();
		dmat subblock(us);
		void setsubblock(us,vd& data); //Set diagonals with this vector
		void setsubblock(us pos,dmat& data); //Set full block
	protected:
		us Nf,Ns,bw;
		dmat Ksub;
		vd Rsub;


};

class eqsystem{
	public:
		eqsystem();
		eqsystem(us Nf,us gp);
		~eqsystem();
		void reset();
		void seteq(us eqnr,us gpnr,subsys& Ksub);
		dmat& getK();
		vd& getR();
	private:
		dmat K; //Sparse matrix
		vd R;
		us Nf,Ns,gp,bw,bh;

};

class equation:public freqoperators::freqoperators{
	public:
		equation();
		equation(us gp,us Nf,d freq);
		subsys Ksub;
	protected:
		us gp,bw;

};

class zeroeq: public equation{
	public:
		zeroeq();
		zeroeq(us gp,us Nf,us varnr);
		subsys getsys(us gpnr);
	protected:
		us varnr;
};

class lincontinuityeq: public equation{
	public:
		lincontinuityeq();
		lincontinuityeq(us gp,us Nf,vd& x,d freq);
		subsys getsys(us gpnr);

		virtual ~lincontinuityeq();
	protected:
		d rho0=1.2;
		d pm=101325;
		vd x;
		d c0sq=pow(343.0,2);
		d K1,K2;
	private:
};
class linmomentumeq: public equation{
	public:
		linmomentumeq();
		linmomentumeq(us gp,us Nf,vd& x,d freq);
		subsys getsys(us gpnr);

		virtual ~linmomentumeq();
	protected:
		d rho0=1.2;
		d pm=101325;
		vd x;
		d c0sq=pow(343.0,2);
		d K1,K2;
	private:
};

class solution{
	public:
		solution();
		solution(us gp,us Nf,vd& sol);
		void setsol(vd& sol);

		vd& getU();
		vd& getp();
		vd getU(us freqnr);
		vd getp(us freqnr);
	protected:
		vd getvarsol(us varnr);
		us Unr=1,rhonr=0,Tnr=2,pnr=3,Tsnr=4;
		us Ns,Nf,gp;
		vd sol;
		vd Usol;
		vd psol;
};

} //namespace eq
#endif // CONTINUITYEQ_H
