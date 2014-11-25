#pragma once
#ifndef _VARNR_H_
#define _VARNR_H_

#define RHONR (0)
#define UNR (1)
#define TNR (2)
#define PNR (3)
#define TSNR (4)
// Number of variables
#define NVARS (5)
namespace tube{

  enum varnr{rho=RHONR,U=UNR,T=TNR,p=PNR,Ts=TSNR};


}                // namespace tube

#endif /* _VARNR_H_ */
