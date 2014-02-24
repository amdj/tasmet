#ifndef HEATEQ_VERTPLATES_H
#define HEATEQ_VERTPLATES_H

#include "var.h"
#include "gas.h"
#include "vtypes.h"
#include "math_anne.h"
using namespace variable;
using namespace gases;

namespace vertplates{




class heateq_vertplates
{
	public:
		heateq_vertplates(d y0,varoperations& fop);
		virtual ~heateq_vertplates();

		void setdata(const var& Tw,const var& T,const var& p,const var& dTdx,const var& dpdx,gas& g);
		vd Tdistr(d t,vd y); //compute the temperature distribution as a function of y-coordinate for a given time

		vd Temp(d time,vd y); //Compute the temperature distribution for a given time and position vector along the y-axis (typically 0<=y<=y0 , or -y0<=y<=y0)
		var Temp(d y); //Compute the temperature for a given y in terms of a var
	private:
		var g_n(d y); //Compute g at a specific y
		vd Cerr(const vd& C);
		d y0,rho0,mu0;
	protected:
		vc hnu_n(d y);

		us Nf,Ns;
		vd s; //Shear wave number for each frequency
		var A,uB,C,rho,mu,Tw_over_T,T;
		const varoperations& vop;
		vc An_over_Cnsq(vc& An,vc& Cn) const; //Problems with computing this number related to singularities.

};
} // namespace heateq_vertplates
#endif // HEATEQ_VERTPLATES_H
