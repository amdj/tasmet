#include "conetube.h"
#include "exception.h"

namespace duct{


  ConeTube::ConeTube(const vd& g,d r1,d r2,bool blapprox) 
    :Geom(g,blapprox),rL(r1),rR(r2){
    TRACE(15,"ConeTube::ConeTube()");
    if(r1 <=0 || r2 <= 0)
      throw MyError("Illegal radius given");
    SL=number_pi*pow(r1,2);
    SR=number_pi*pow(r2,2);
    setPrismatic(false);
  }
  d ConeTube::S(us i) const{
    return SL+(SR-SL)*x(i)/L();
  }
  d ConeTube::rh(us i) const{
    return 0.5*(rL+(rR-rL)*x(i)/L());
  }
  void ConeTube::show() const{
    cout << "Conical tube\n";
    cout << "Number of cells:  " << nCells() << "\n";
    cout << "Length:           " << L() << "\n";
    cout << "Radius left side: " << rL<<"\n";
    cout << "Radius right side:" << rR<<"\n";
  }


  CylindricalTube::CylindricalTube(const vd& g,d r,bool blapprox):ConeTube(g,r,r,blapprox){
    setPrismatic(true);
  }
  void CylindricalTube::show() const{
    cout << "Cylindrical tube\n";
    cout << "Number of cells: " << nCells() << "\n";
    cout << "Length:          " << L() << "\n";
    cout << "Radius:          " << rL << "\n";
  }



} // namespace duct


