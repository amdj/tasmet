#include "globalconf.h"


namespace tasystem{
  
  Globalconf::Globalconf(us Nf,d freq): Globalconf(Nf,freq,"air"){}
  Globalconf::Globalconf(us Nf,d freq,string Gas):Globalconf(Nf,freq,Gas,293.15,101325,0){}
  Globalconf::Globalconf(us Nf,d freq,string gas,d T0,d p0,d Mass):
    Nf(Nf),
    Ns(2*Nf+1),
    gas(gas),
    freq(freq),
    omg(2*pi*freq),
    T0(T0),
    p0(p0),
    Mass(Mass)
{
  TRACE(10,"Globalconf constructor done");
}
  Globalconf::~Globalconf(){}


} // Namespace tasystem
