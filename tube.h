/*
 * lintube.h
 *
 *  Created on: Oct 8, 2013
 *      Author: anne
 */

#ifndef TUBE_H_
#define TUBE_H_
#include "vvar.h"
#include "vtypes.h"
#include <string.h>
namespace tube {

class tube {
	protected:
		dmat Ddt,fF,iF;

	public:
		us Ns,Nf,gp;
		vvar::vvar rho,p,U,u,T,rhoU,rhoUu;//,Etot,Hf;
		us DofsPerVar;
		tube(us,us);// gp,Nf
		virtual ~tube();
};

} /* namespace tube */
#endif /* TUBE_H_ */
