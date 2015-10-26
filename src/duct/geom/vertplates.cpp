#define PERCENT (1.0/100)
// #include "geomhelpers.h"
#include "vertplates.h"
#include "vtypes.h"
#include "bessel.h"
#include "skewsine.h"
#include <cassert>

namespace duct{


  VertPlates::VertPlates(const vd& g,d S,d phi,d y0,bool blapprox):
    Geom(g,blapprox,true),
    S_(S),
    phi_(phi),
    rh_(y0)
  {
    TRACE(15,"VertPlates::VertPlates()");
    assert(y0>0);
    assert(0<phi && phi<=1.0);
  }
  void VertPlates::show() const{
    cout << "VertPlates geom\n";
    cout << "Number of cells:   " << nCells() << "\n";
    cout << "Length:            " << L() << "\n";
    cout << "Hydraulic radius:  " << rh_<<"\n";
    cout << "Total surface area:" << S_<<"\n";
    cout << "Porosity:          " << phi_<<"\n";
  }

} // namespace duct

