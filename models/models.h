#ifndef _MODELS_H_
#define _MODELS_H_

#include "solver.h"
#include "tasystem.h"
using tasystem::Solver;
SPOILNAMESPACE
// Solver* ConeTube(us gp,us Nf,d freq,d L,d r1,d r2,vd p1,int loglevel,d kappa,us isentropic,us blapprox);

#define ISENTROPIC 1
#define BLAPPROX 2
#define DRIVEN 4
#define ISOTWALL 8



Solver* SimpleTube(us gp,us Nf,d freq,d L,d r,d Tl,d Tr,vd p1,int loglevel,d kappa,int options);
Solver* Fubini(us gp,us Nf,d freq,d L,vd p1,int loglevel,d kappa,int options);
Solver* ThreeTubes(us gp,us Nf,d freq,d L,d R1,d R2,vd p1,int loglevel,d kappa,d Tr,int options);
Solver* Atchley_Engine(us gp,us Nf,d freq,d Tr,int loglevel,d kappa,vd p1,d p0,int options);

#endif /* _MODELS_H_ */

