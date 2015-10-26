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
  void Globalconf::set(us Nf,d freq){
    TRACE(15,"Globalconf::set(Nf_,freq)");
    //ctor
    this->Nf_=Nf;
    this->omg=2*number_pi*freq;
    us Ns=this->Ns();
    // Reinitialize all operators
    iDFT_=zeros<dmat>(Ns,Ns);
    fDFT_=zeros<dmat>(Ns,Ns);

    DDTfd_=zeros<dmat>(Ns,Ns);
    fDFT_.row(0).fill(1.0/double(Ns));

    for(us i=1;i<=Nf_;i++){
      for(us j=0; j<Ns;j++){
        //Row i+1 (cosine components)
        fDFT_(2*i-1,j)=2.0*cos(2.0*number_pi*double(i)*double(j)/double(Ns))/double(Ns);
        //Row i (sine components)
        fDFT_(2*i,j)=-2.0*sin(2.0*number_pi*double(i)*double(j)/double(Ns))/double(Ns);
      }
    }
    iDFT_.col(0).fill(1.0);	// Steady part
    for(us k=0;k<Ns;k++){
      for (us n=1;n<=Nf_;n++){
        iDFT_(k,2*n-1)=cos(2.0*number_pi*double(n)*double(k)/Ns);
        iDFT_(k,2*n)=-sin(2.0*number_pi*double(n)*double(k)/Ns);
      }
    }
    for(us i=1;i<=Nf_;i++){
      DDTfd_(2*i-1,2*i  )=-double(i)*omg;
      DDTfd_(2*i  ,2*i-1)=double(i)*omg;
    }

  }
  void Globalconf::setNf(us Nf){set(Nf,getfreq());}  
  void Globalconf::setfreq(d freq){set(Nf(),freq);}
  void Globalconf::setomg(d omg)  {set(Nf(),omg/(2*number_pi));}

} // Namespace tasystem
