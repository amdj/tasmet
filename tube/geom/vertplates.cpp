#define PERCENT (1.0/100)
// #include "geomhelpers.h"
#include "vertplates.h"
#include "vtypes.h"
#include "bessel.h"
#include "skewsine.h"

namespace tube{

  VertPlates::VertPlates(const Grid& g,d S,d phi,d y0,bool blapprox):
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
  void TransitionVertPlates::show() const{
    cout << "TransitionVertPlates geom\n";
    cout << "Number of cells:   " << nCells() << "\n";
    cout << "Length:            " << L() << "\n";
    cout << "Hydraulic radius:  " << rh_<<"\n";
    cout << "Total surface area:" << S_<<"\n";
    cout << "Porosity:          " << phi_<<"\n";
    cout << "Transition surface area: " << S_other << "\n";
    cout << "This Surface area      : " << S_ << "\n";
    if(transition.Position()==left)
      cout << "Smooth transition made on the left side\n";
    else
      cout << "Smooth transition made on the right side\n";
    cout << "Percentage of tube length involved with transition: " << perc<<"\n";
  }

  Transition::Transition(pos position,d perc):position(position),perc(perc){}
  d Transition::percd_other(d x_ov_L) const {
    using math_common::skewsine;
    d percd=perc/100;           // percentage as decimal
    if(position==left){
      if(x_ov_L<percd)
        return 1-skewsine(x_ov_L/percd);
      else
        return 0;
    } // position==left
    else{
      if(x_ov_L>1-percd)
        return skewsine((x_ov_L-(1-percd))/percd);
      else
        return 0;
    } // position==right

  }
  TransitionVertPlates::TransitionVertPlates(const Grid& grid,\
                                             d S,d phi,d y0,
                                             pos TransitionSide,    \
                                             const Geom& other,\
                                             pos sideofremote,\
                                             d perc,\
                                             bool blapprox):
    VertPlates(grid,S,phi,y0,blapprox),
    transition(TransitionSide,perc)
  {
    TRACE(14,"TransitionVertPlates::TransitionVertPlates()");
    setPrismatic(false);
    if(sideofremote==left)    {
      S_other=other.Sleft();
      phi_other=other.phileft();
      WARN("rh interpolated as well!");
      rh_other=other.rhleft();
    }
    else {
      S_other=other.Sright();
      phi_other=other.phiright();
      WARN("rh interpolated as well!");
      rh_other=other.rhright();
    }
  }

  d TransitionVertPlates::S(us i) const{
    TRACE(5,"TransitionVertPlates::S()");
    d xi=x(i);
    d L=this->L();
    using math_common::skewsine;
    assert(i<gp());
    d percd=transition.percd_other(xi/L);
    return (1-percd)*S_+percd*S_other;
  }
  d TransitionVertPlates::phi(us i) const{
    TRACE(5,"TransitionVertPlates::phi()");
    d xi=x(i);
    d L=this->L();
    using math_common::skewsine;
    assert(i<gp());
    d percd=transition.percd_other(xi/L);
    return (1-percd)*phi_+percd*phi_other;
  }
  d TransitionVertPlates::rh(us i) const{
    TRACE(5,"TransitionVertPlates::rh()");
    d xi=x(i);
    d L=this->L();
    using math_common::skewsine;
    assert(i<gp());
    d percd=transition.percd_other(xi/L);
    return (1-percd)*rh_+percd*rh_other;
  }


} // namespace tube




