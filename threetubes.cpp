#include "models.h"
SPOILNAMESPACE
#include "tube.h"

int main(int argc,char* argv[]) {
  cout <<  "Running test..." << endl;
  int loglevel=20;
  us Nf=0;
  us gp=4;
  d freq=500;
  
  if(argc>1)
    loglevel=atoi(argv[1]);
  if(argc>2)
    Nf=atoi(argv[2]);
  if(argc>3)
    gp=atoi(argv[3]);
  cout << "gp:"<< gp<<"\n";
  // us gp=40;
  // us Nf=1;
  // us Ns=2*Nf+1;
  // double f=100;
  // double omg=2*number_pi*f;
  // double T=1/f;
  d T0=293.15;
  cout<< "Loglevel:"<<loglevel<<"\n";
  initlog(loglevel);
  vd p1(2*Nf+1,arma::fill::zeros);
  if(Nf>0)
    p1(1)=0.5e-2;
  d r=0.5e-2;
  d S1=number_pi*pow(r,2);
  d S2=S1;
  d kappa=1;
  d L=0.1;
  Solver* sol=ThreeTubesEngine(gp,Nf,freq,T0+10,p1,loglevel,kappa);

  // sol->sys().show(true);
  // sol->doIter();
  // sol->doIter();
  // cout <<sol->sys().jac();
  sol->solve();
  // cout <<"Result:\n"<<  sol->sys.getRes();
  cout <<"Result:\n"<< static_cast<tube::Tube*>(sol->sys().getSeg(0))->getResAt(2,0);
  cout <<"Result:\n"<< static_cast<tube::Tube*>(sol->sys().getSeg(1))->getResAt(2,0);
  cout <<"Result:\n"<< static_cast<tube::Tube*>(sol->sys().getSeg(2))->getResAt(2,0);

  
  
  
  delete sol;
  return 0;
}








