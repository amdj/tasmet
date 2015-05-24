#include "globalconf.h"
#include "constants.h"
#include "exception.h"

namespace tasystem{
  Globalconf Globalconf::airSTP(us Nf_,d freq) {
    return Globalconf(Nf_,freq,"air",constants::T0,constants::p0);
  }
  Globalconf Globalconf::heliumSTP(us Nf_,d freq) {
    return Globalconf(Nf_,freq,"helium",constants::T0,constants::p0);
  }
  Globalconf::Globalconf(us Nf,d freq,const string& gasstring,d T0,d p0):
    gas_(gasstring)
  {
    // Sanity checks
    if(2*number_pi*freq<constants::minomg && 2*number_pi*freq>constants::maxomg)
      throw MyError("Illegal frequency given");
    if(Nf>=constants::maxNf)
      throw("Too large number of frequencies given");
    if(T0>constants::maxT || T0< constants::minT)
      throw MyError("Illegal reference temperature given");
    if(p0>constants::maxp || p0< constants::minp)
      throw MyError("Illegal reference pressure given");
    // End sanity checks

    p0_=p0;
    T0_=T0;
    // Initialize FFT matrices
    set(Nf,freq);
    TRACE(10,"Globalconf constructor done");
  }
  
  void Globalconf::show() const {
    cout << "------- Globalconf configuration ------ \n";        
    cout << "------- Nf             : "<< Nf_ <<"\n";
    cout << "------- Base frequency : " << omg/2/number_pi << " Hz\n";           
    cout << "------- Gas            : " << string(gas_) << "\n";          
    cout << "------- p0             : " << p0_ << " [Pa] \n";                 
    cout << "------- T0             : " << T0_ << " [K] \n";
    cout << "------- rho0           : " << rho0() << " [kg/m^3] \n";
    cout << "------- c0:            : " << c0() << "\n";

  }
  void Globalconf::setNf(us Nf){
    set(Nf,this->omg/2/number_pi);
  }
  void Globalconf::set(us Nf,d freq){
    TRACE(15,"Globalconf::set(Nf_,freq)");
    //ctor
    this->Nf_=Nf;
    Ns_=2*Nf_+1;
    // Reinitialize all operators
    iDFT=zeros<dmat>(Ns_,Ns_);
    fDFT=zeros<dmat>(Ns_,Ns_);

    DDTfd=zeros<dmat>(Ns_,Ns_);
    DDTtd=zeros<dmat>(Ns_,Ns_);
    ddt=zeros<dmat>(Ns_-1,Ns_-1);
    iddt=zeros<dmat>(Ns_-1,Ns_-1);
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
    fDFT.row(0).fill(1.0/double(Ns_));

    for(us i=1;i<=Nf_;i++){
      for(us j=0; j<Ns_;j++){
        //Row i+1 (cosine components)
        fDFT(2*i-1,j)=2.0*cos(2.0*number_pi*double(i)*double(j)/double(Ns_))/double(Ns_);
        //Row i (sine components)
        fDFT(2*i,j)=-2.0*sin(2.0*number_pi*double(i)*double(j)/double(Ns_))/double(Ns_);
      }
    }
  }
  void Globalconf::updateiDFT(){
    TRACE(10,"Globalconf::updateiDFT()");
    
    iDFT.col(0).fill(1.0);	// Steady part
    for(us k=0;k<Ns_;k++){
      for (us n=1;n<=Nf_;n++){
        iDFT(k,2*n-1)=cos(2.0*number_pi*double(n)*double(k)/Ns_);
        iDFT(k,2*n)=-sin(2.0*number_pi*double(n)*double(k)/Ns_);
      }
    }
  }
  void Globalconf::updateiomg(){
    TRACE(15,"Globalconf::updateiomg()");
    if(Nf_!=0){
      for(us i=1;i<=Nf_;i++){
        DDTfd(2*i-1,2*i  )=-double(i)*omg;
        DDTfd(2*i  ,2*i-1)=double(i)*omg;
      }
      DDTtd=iDFT*DDTfd*fDFT;
      ddt=DDTfd.submat(1,1,Ns_-1,Ns_-1);
      //cout << "ddt:" << ddt << endl;
      iddt=inv(ddt);
    }
    else{
      DDTfd(0,0)=0;
    }
  }
  



} // Namespace tasystem
