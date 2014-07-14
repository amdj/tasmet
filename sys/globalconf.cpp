#include "globalconf.h"

LOGTHIS
namespace tasystem{
  
  Globalconf::Globalconf(us Nf,d freq,string gas,d T0,d p0,d Mach,d S0,d dx,d Mass,d kappa):
    gas(gas)
  {
    this->Gastype=gas;
    this->p0=p0;
    this->T0=T0;
    this->S0=S0;
    this->dx=dx;
    this->V0=S0*dx;
    this->Mass=Mass;
    this->c0=this->gas.cm(T0);
    this->kappa=kappa;
    this->Mach=Mach;
    if(Nf==0 || Mach<1e-10)
      M=1.0;
    else
      this->M=Mach;
    set(Nf,freq);
    TRACE(10,"Globalconf constructor done");
    
  }
  void Globalconf::show(){
    cout << "------- Globalconf configuration ------ \n"			\
	 << "------- Nf             : "<< Nf <<"\n"				\
	 << "------- Base frequency : " << freq << " Hz\n"			\
	 << "------- Gas            : " << Gastype << "\n"			\
	 << "------- p0             : " << p0 << " [Pa] \n"			\
	 << "------- T0             : " << T0 << " [K] \n"			\
	 << "------- kappa:         : " << kappa << "\n"			\
      ;
    

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
    setfreq(freq);

    TRACE(-1,"fDFT:" << fDFT);
  }
  void Globalconf::setfreq(d freq)  {
    this->freq=freq;
    oldomg=omg;
    omg=2.0*number_pi*freq;
    updateiDFT();
    updatefDFT();
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
