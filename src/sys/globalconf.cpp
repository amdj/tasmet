#include "globalconf.h"

TRACETHIS
namespace tasystem{
  Globalconf Globalconf::airSTP(us Nf__,d freq,d kappa){
    return Globalconf(Nf__,freq,"air",293.15,101325,kappa);
  }
  Globalconf::Globalconf(us Nf__,d freq,const string& gasstring,d T0,d p0    \
                         ,d kappa,bool driven):
    driven(driven),
    gas(gasstring)
  {
    // Sanity checks
    assert(2*number_pi*freq>MINOMG && 2*number_pi*freq<MAXOMG);
    assert(Nf__<MAXNF);
    assert(T0<2000 && T0>0);
    assert(p0>0);
    assert(Mass>=0);
    assert(kappa>0);
    // End sanity checks

    this->p0=p0;
    this->T0=T0;
    this->kappa=kappa;

    // Initialize FFT matri
    set(Nf__,freq);
    TRACE(10,"Globalconf constructor done");
    
  }
  
  void Globalconf::show() const {
    cout << "------- Globalconf configuration ------ \n"        \
         << "------- Nf             : "<< Nf__ <<"\n"                       \
         << "------- Base frequency : " << omg/2/number_pi << " Hz\n"       \
         << "------- Gas            : " << getGas() << "\n"                 \
         << "------- p0             : " << p0 << " [Pa] \n"                 \
         << "------- T0             : " << T0 << " [K] \n"                  \
         << "------- rho0           : " << rho0() << " [kg/m^3] \n"         \
         << "------- kappa:         : " << kappa << "\n"                    \
         << "------- c0:            : " << c0() << "\n" ;

  }
  void Globalconf::setNf(us Nf){
    set(Nf,this->omg/2/number_pi);
  }
  void Globalconf::set(us Nf__,d freq){
    TRACE(15,"Globalconf::set(Nf__,freq)");
    //ctor
    this->Nf__=Nf__;
    Ns__=2*Nf__+1;
    // Reinitialize all operators
    iDFT=zeros<dmat>(Ns__,Ns__);
    fDFT=zeros<dmat>(Ns__,Ns__);

    DDTfd=zeros<dmat>(Ns__,Ns__);
    DDTtd=zeros<dmat>(Ns__,Ns__);
    ddt=zeros<dmat>(Ns__-1,Ns__-1);
    iddt=zeros<dmat>(Ns__-1,Ns__-1);
    updateiDFT();
    updatefDFT();

    setfreq(freq);

    TRACE(-1,"fDFT:" << fDFT);
  }
  void Globalconf::setfreq(d freq){setomg(2*number_pi*freq);}
  void Globalconf::setomg(d omg)  {
    TRACE(15,"Globalconf::setomg()");
    
    this->omg=omg;
    oldomg=omg;
    updateiomg();
  }
  void Globalconf::updatefDFT(){
    fDFT.row(0).fill(1.0/double(Ns__));

    for(us i=1;i<=Nf__;i++){
      for(us j=0; j<Ns__;j++){
        //Row i+1 (cosine components)
        fDFT(2*i-1,j)=2.0*cos(2.0*number_pi*double(i)*double(j)/double(Ns__))/double(Ns__);
        //Row i (sine components)
        fDFT(2*i,j)=-2.0*sin(2.0*number_pi*double(i)*double(j)/double(Ns__))/double(Ns__);
      }
    }
  }
  void Globalconf::updateiDFT(){
    TRACE(10,"Globalconf::updateiDFT()");
    
    iDFT.col(0).fill(1.0);	// Steady part
    for(us k=0;k<Ns__;k++){
      for (us n=1;n<=Nf__;n++){
        iDFT(k,2*n-1)=cos(2.0*number_pi*double(n)*double(k)/Ns__);
        iDFT(k,2*n)=-sin(2.0*number_pi*double(n)*double(k)/Ns__);
      }
    }
  }
  void Globalconf::updateiomg(){
    TRACE(15,"Globalconf::updateiomg()");
    if(Nf__!=0){
      for(us i=1;i<=Nf__;i++){
        DDTfd(2*i-1,2*i  )=-double(i)*omg;
        DDTfd(2*i  ,2*i-1)=double(i)*omg;
      }
      DDTtd=iDFT*DDTfd*fDFT;
      ddt=DDTfd.submat(1,1,Ns__-1,Ns__-1);
      //cout << "ddt:" << ddt << endl;
      iddt=inv(ddt);
    }
    else{
      DDTfd(0,0)=0;
    }
    omgvec=vd(Nf__+1);
    for(us i=0; i<Nf__+1;i++)
      omgvec(i)=omg*i;
  }
  



} // Namespace tasystem
