#define PERCENT (1.0/100)
// #include "geomhelpers.h"
#include "vertplates.h"
#include "vtypes.h"
#include "bessel.h"

namespace geom{

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


  // void smoothEnds(Geom& smooththis,int pos,
  //                 const Geom& to,int topos,int perc){
  //   TRACE(3,"SmoothEnd()");
  //   d dx=PERCENT*perc*smooththis.L; // Length of segment
  //   // to smooth out
  //   TRACE(15,"Perc: "<< perc);
  //   TRACE(15,"First L: "<< smooththis.L);
  //   // TRACE(15,"Second L: "<< second.L<<"\n" );        
  //   TRACE(15,"dx to smooth out:"<< dx);

  //   Geom tempgeom=smooththis;    
  //   // string fcshape=first.shape;

    
  //   us i,j;    
  //   if(pos==FIRST)      {
  //     i=0;
  //     TRACE(15,"Smoothing beginning of first Geom, i="<< i);
  //   }
  //   else{
  //     i=tempgeom.nCells;
  //     TRACE(15,"To end of first Geom, i="<< i);
  //   }
  //   if(topos==FIRST){
  //     j=0; 
  //     TRACE(15,"Smoothing beginning of second Geom, j="<< j);
  //   }
  //   else{
  //     j=to.nCells;
  //     TRACE(15,"To end of second Geom, j=" << j);
  //   }

  //   // Now adjust it
  //   TRACE(15,"SFSG");
  //   d Sorig=smooththis.S(i);
  //   d Sto=to.S(j);

  //   d phiorig=smooththis.phi(i);
  //   d phito=to.phi(j);

  //   d rhorig=smooththis.rh(i);
  //   d rhto=to.rh(j);

  //   TRACE(15,"SFSG");
  //   if(pos==FIRST)      {
  //     TRACE(15,"SFSG");
          
  //     us i=0;
  //     while(tempgeom.x(i)<dx){
  //       tempgeom.S(i)=Sto+(Sorig-Sto)*math_common::skewsine(tempgeom.x(i)/dx);
  //       tempgeom.phi(i)=phito+(phiorig-phito)*math_common::skewsine(tempgeom.x(i)/dx);
  //       tempgeom.rh(i)=rhto+(rhorig-rhto)*math_common::skewsine(tempgeom.x(i)/dx);          
  //       i++;
  //     }
  //   }
  //   else{
  //     TRACE(15,"SFSG");
  //     us ilast=smooththis.x.size()-1;
  //     d L=smooththis.x(ilast);
  //     us i=ilast;
  //     TRACE(15,"ilast:"<< ilast)
  //     while(tempgeom.x(i)>(L-dx)){
  //       TRACE(15,"i="<<i);
  //       tempgeom.S(i)=Sto+(Sorig-Sto)*math_common::skewsine((L-smooththis.x(i))/dx);
  //       tempgeom.phi(i)=phito+(phiorig-phito)*math_common::skewsine((L-smooththis.x(i))/dx);
  //       tempgeom.rh(i)=rhto+(rhorig-rhto)*math_common::skewsine((L-smooththis.x(i))/dx);
  //       i--;
  //     }
  //   }
  //   smooththis=tempgeom;
  // }
} // namespace geom




