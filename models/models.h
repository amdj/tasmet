#ifndef _MODELS_H_
#define _MODELS_H_

#include "solver.h"
using tasystem::Solver;
SPOILNAMESPACE
// Solver* ConeTube(us gp,us Nf,d freq,d L,d r1,d r2,vd p1,int loglevel,d kappa,us isentropic,us blapprox);


Solver* SimpleTube(us gp,us Nf,d freq,d L,d r,d Tl,d Tr,vd p1,int loglevel,d kappa,us blapprox=1,us driven=1,us rwall=0);
Solver* Fubini(us gp,us Nf,d freq,d L,vd p1,int loglevel,d kappa,us fullenergy=0);
Solver* ThreeTubes(us gp,us Nf,d freq,d L,d R1,d R2,vd p1,int loglevel,d kappa,d Tr,us isentropic);
Solver* Atchley_Engine(us gp,us Nf,d freq,d Tr,int loglevel,d kappa,us driven,vd p1,d p0);

#endif /* _MODELS_H_ */

