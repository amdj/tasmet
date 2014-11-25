#include "conetube.h"


namespace geom{

  ConeTube::ConeTube(const Grid& g,d r1,d r2,bool blapprox):Geom(g,blapprox),rL(r1),rR(r2){
    TRACE(15,"ConeTube::ConeTube()");
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


  CylindricalTube::CylindricalTube(const Grid& g,d r,bool blapprox):ConeTube(g,r,r,blapprox){
    setPrismatic(true);
  }
  void CylindricalTube::show() const{
    cout << "Cylindrical tube\n";
    cout << "Number of cells: " << nCells() << "\n";
    cout << "Length:          " << L() << "\n";
    cout << "Radius:          " << rL << "\n";
  }

  TransitionCylindricalTube::TransitionCylindricalTube(const Grid& g,d r,\
                                                       pos transitionside,\
                                                       const Geom& other, \
                                                       pos connectionsideother,\
                                                       d perc,bool blapprox):
    Geom(g,blapprox),
    transition(transitionside,perc),
    r_(r)
  {
    TRACE(15,"TransitionCylindricalTube::TransitionCylindricalTube()");

    setPrismatic(false);
    S_=number_pi*pow(r,2);
    if(connectionsideother==left)    {
      S_other=other.Sleft();
      phi_other=other.phileft();
    }
    else {
      S_other=other.Sright();
      phi_other=other.phiright();
    }
  } // Constructor
  d TransitionCylindricalTube::S(us i) const{
    TRACE(5,"TransitionCylindricalTube::S()");
    d xi=x(i);
    d L=this->L();
    assert(i<gp());
    d percd=transition.percd_other(xi/L);
    return (1-percd)*S_+percd*S_other;
  }

  d TransitionCylindricalTube::phi(us i) const {
    TRACE(5,"TransitionCylindricalTube::phi()");
    d xi=x(i);
    d L=this->L();
    d phi_=1.0;
    assert(i<gp());
    d percd=transition.percd_other(xi/L);
    return (1-percd)*phi_+percd*phi_other;
  }
  void TransitionCylindricalTube::show() const{
    cout << "Type: TransitionCylindricalTube\n";
    cout << "Number of cells        : " << nCells() << "\n";
    cout << "Length                 : " << L() << "\n";
    cout << "Radius of tube to go to: " << r_<<"\n";
    cout << "This Surface area      : " << S_ << "\n";
    cout << "Transition surface area: " << S_other << "\n";
    cout << "\n";
    if(transition.Position()==left)
      cout << "Smooth transition made on the left side\n";
    else
      cout << "Smooth transition made on the right side\n";
    cout << "Percentage of tube length involved with transition: " << perc<<"\n";
  }



} // namespace geom


