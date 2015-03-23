#pragma once
#ifndef _VARNR_H_
#define _VARNR_H_

#ifndef SWIG
const int RHONR=0;
const int UNR=1;
const int TNR=2;
const int PNR=3;
const int TSNR=4;
// Number of variables
const int NVARS=5;
#endif


namespace tube{

  enum varnr{rho=RHONR,
             U=UNR,
             T=TNR,
             p=PNR,
             Ts=TSNR,
             rhoL,
             rhoR,
             UL,
             UR,
             TL,
             TR,
             pL,
             pR,
             TsL,
             TsR
  };

  enum physquant{massFlow,
                 momentumFlow,
                 energyFlow,
                 heatFlow,
                 solidHeatFlow,
                 rhoRT,
  };
}                // namespace tube

#endif /* _VARNR_H_ */
