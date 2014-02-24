#ifndef LINSOLVER_H
#define LINSOLVER_H
#include "gd.h"
#include "lintube.h"
#include "lintubefv.h"

namespace solver{
using std::cout;
using std::endl;
class linsolver
{

	protected:


	public:

		vd getResult();
		vd getpResult();
		vd getx();
		vd getV();
		vd getUResult();
		void setResult(vd);
		vd getError();
		vd getUError();
		vd getUResult(int);
		vd getpResult(int);
		vd getpError();
		virtual ~linsolver();
		void setSubVec(vd*,int,vd);
		vd getSubVec(vd bigvec,int s,int e);
		void setpEr(double*);
		//boost::python::object stdVecToNumpyArray(vd const&);
};
}
#endif // LINSOLVER_H
