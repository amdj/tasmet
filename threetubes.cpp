#include "models.h"
SPOILNAMESPACE

int main(int argc,char* argv[]) {
  cout <<  "Running test..." << endl;
  int loglevel=25;
  us  Nf=4;
  us gp=60;
  d freq=100;
  
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
  Solver* sol=ThreeTubes(gp,Nf,freq,L,S1,S2,p1,loglevel,kappa);

  sol->sys.show(true);

  sol->solve();
  // sol->sys.getRes();
  // sol->sys->Init();
  // cout << "segjac\n" << sol->sys.getSeg(0)->Jac();
  // cout << "\n" << sol->sys->Jac();
  // dmat jac=sol->sys.jac();
  // cout << "error\n" << sol->sys.error();
  // cout << "Jac:\n" <<jac;

  // cout << "\nDeterminant jac:"<< arma::det(jac)<<"\n"; 
  // sol->doIter();
  // sol->sys->show(true);

  // cout << sol->sys.getRes();
  delete sol;
  return 0;
}








