#include "globalconf.h"

TRACETHIS
namespace tasystem{
  Globalconf Globalconf::airSTP(us Nf,d freq,d kappa){
    return Globalconf(Nf,freq,"air",293.15,101325,kappa);
  }
  Globalconf::Globalconf(us Nf,d freq,string gas,d T0,d p0,d kappa):
    gas(gas)
  {
    // Sanity checks
    assert(2*number_pi*freq>MINOMG && 2*number_pi*freq<MAXOMG);
    assert(Nf<MAXNF);
    assert(T0<2000 && T0>0);
    assert(p0>0);
    assert(Mass>=0);
    assert(kappa>0);
    // End sanity checks

    this->Gastype=gas;
    this->p0=p0;
    this->T0=T0;
    this->rho0=this->gas.rho(T0,p0);

    this->c0=this->gas.cm(T0);
    this->kappa=kappa;

    // Initialize FFT matri
    set(Nf,freq);
    TRACE(10,"Globalconf constructor done");
    
  }
  void Globalconf::show() const {
    cout << "------- Globalconf configuration ------ \n"			\
	 << "------- Nf             : "<< Nf <<"\n"				\
	 << "------- Base frequency : " << omg/2/number_pi << " Hz\n"			\
	 << "------- Gas            : " << Gastype << "\n"			\
	 << "------- p0             : " << p0 << " [Pa] \n"			\
	 << "------- T0             : " << T0 << " [K] \n"			\
	 << "------- rho0           : " << rho0 << " [kg/m^3] \n"		\
	 << "------- kappa:         : " << kappa << "\n"			\
	 << "------- c0:            : " << c0 << "\n" ;
    

  }
  void Globalconf::set(us Nf,d freq){
    TRACE(0,"Globalconf::set(Nf,freq)");
    //ctor
    this->Nf=Nf;
    Ns=2*Nf+1;
    // Reinitialize all operators
    iDFT=zeros<dmat>(Ns,Ns);
    fDFT=zeros<dmat>(Ns,Ns);

    DDTfd=zeros<dmat>(Ns,Ns);
    DDTtd=zeros<dmat>(Ns,Ns);
    ddt=zeros<dmat>(Ns-1,Ns-1);
    iddt=zeros<dmat>(Ns-1,Ns-1);
    updateiDFT();
    updatefDFT();

    setfreq(freq);

    TRACE(-1,"fDFT:" << fDFT);
  }
  void Globalconf::setfreq(d freq){setomg(2*number_pi*freq);}
  void Globalconf::setomg(d omg)  {
    this->omg=omg;
    oldomg=omg;
    updateiomg();
  }
  void Globalconf::updatefDFT(){
    fDFT.row(0).fill(1.0/double(Ns));

    for(us i=1;i<=Nf;i++){
      for(us j=0; j<Ns;j++){
	//Row i+1 (cosine components)
	fDFT(2*i-1,j)=2.0*cos(2.0*number_pi*double(i)*double(j)/double(Ns))/double(Ns);
	//Row i (sine components)
	fDFT(2*i,j)=-2.0*sin(2.0*number_pi*double(i)*double(j)/double(Ns))/double(Ns);
      }
    }
  }
  void Globalconf::updateiDFT(){
    iDFT.col(0).fill(1.0);	// Steady part
    for(us k=0;k<Ns;k++){
      for (us n=1;n<=Nf;n++){
	iDFT(k,2*n-1)=cos(2.0*number_pi*double(n)*double(k)/Ns);
	iDFT(k,2*n)=-sin(2.0*number_pi*double(n)*double(k)/Ns);
      }
    }
  }
  void Globalconf::updateiomg(){
    TRACE(0,"Globalconf::updateiomg()");
    if(Nf!=0){
      for(us i=1;i<=Nf;i++){
	DDTfd(2*i-1,2*i  )=-double(i)*omg;
	DDTfd(2*i  ,2*i-1)=double(i)*omg;
      }
      DDTtd=iDFT*DDTfd*fDFT;
      ddt=DDTfd.submat(1,1,Ns-1,Ns-1);
      //cout << "ddt:" << ddt << endl;
      iddt=inv(ddt);
    }
    else{
      DDTfd(0,0)=0;
    }
    omgvec=vd(Nf+1);
    for(us i=0; i<Nf+1;i++)
      omgvec(i)=omg*i;
  }
  



} // Namespace tasystem