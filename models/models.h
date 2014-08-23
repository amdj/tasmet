#ifndef _MODELS_H_
#define _MODELS_H_

#include "solver.h"
using tasystem::Solver;
SPOILNAMESPACE
Solver* ConeTube(us gp,us Nf,d freq,d L,d r1,d r2,vd p1,int loglevel,d kappa);
Solver* ThreeTubes(us gp,us Nf,d freq,d L,d S1,d S2,vd p1,int loglevel,d kappa);
Solver* Fubini(us gp,us Nf,d freq,d L,d S,vd p1,int loglevel,d kappa);
Solver* Fubini_fullenergy(us gp,us Nf,d freq,d L,d S,vd p1,int loglevel,d kappa);
Solver* ThreeTubesConduction(us gp,us Nf,d freq,d L,d S1,d S2,vd p1,int loglevel,d kappa,d Tr);

Solver* ThreeTubesEngineDriven(us gp,us Nf,d freq,d Tr,vd p1,int loglevel,d kappa);
Solver* ThreeTubesEngine(us gp,us Nf,d freq,d Tr,int loglevel,d kappa);

#endif /* _MODELS_H_ */

