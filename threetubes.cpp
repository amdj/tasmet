#include "models.h"
SPOILNAMESPACE
#include "tube.h"

int main(int argc,char* argv[]) {
  cout <<  "Running test..." << endl;
  int loglevel=20;
  us Nf=0;
  us gp=4;
  d freq=500;
  d p0=101325;
  if(argc>1)
    loglevel=atoi(argv[1]);
  if(argc>2)
    Nf=atoi(argv[2]);
  if(argc>3)
    gp=atoi(argv[3]);
  inittrace(loglevel);
  cout << "gp:"<< gp<<"\n";

  d T0=293.15;
  d Tr=T0;
  int options=0;
  // options|=DRIVEN|ISENTROPIC;
  cout<< "Loglevel:"<<loglevel<<"\n";

  vd p1(2*Nf+1,arma::fill::zeros);
  if(Nf>0)
    p1(1)=0.5e-2;
  d r=0.5e-2;
  d R1=r;
  d R2=r;
  d S1=number_pi*pow(r,2);
  d S2=S1;
  d kappa=1;
  d L=0.1;
  Solver* sol=ThreeTubes(gp,Nf,freq,p0,L,R1,R2,p1,loglevel,kappa,Tr,options);
  sol->sys().show(false);
  // sol->doIter();
  // cout <<sol->sys().jac();
  sol->sys().showJac();
  sol->solve();
  // cout <<"Result:\n"<<  sol->sys.getRes();
  // cout <<"Result:\n"<< static_cast<tube::Tube*>(sol->sys().getSeg(0))->getResAt(2,0);
  // cout <<"Result:\n"<< static_cast<tube::Tube*>(sol->sys().getSeg(1))->getResAt(2,0);
  // cout <<"Result:\n"<< static_cast<tube::Tube*>(sol->sys().getSeg(2))->getResAt(2,0);

  
  
  
  delete sol;
  return 0;
}








