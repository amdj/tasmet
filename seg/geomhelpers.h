#pragma once

#ifndef _GEOMHELPERS_H_
#define _GEOMHELPERS_H_

inline double min(double x,double y)  {
  return x<=y? x : y;
}
inline double max(double x,double y)  {
  return x<=y? y : x;
}


namespace segment{
  const int FIRST=0;
  const int LAST=1;
  class Geom;
  // Smooth geometries of segments together. The percentage is a
  // length part of the smallest of the two given geometries.
  void smoothEnds(Geom& smooththisone,int smooththisonepos,
                  const Geom& tothisone,int tothisonepos,int perc);
  
} // namespace segment

#endif /* _GEOMHELPERS_H_ */
