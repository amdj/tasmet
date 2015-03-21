%module(docstring="TMTubes nonlinear thermoacoustic code") nonlin
%{
  #define PY_ARRAY_UNIQUE_SYMBOL npy_array
  #define SWIG_FILE_WITH_INIT
  #include "pos.h"
  #include "var.h"
  #include "tasystem.h"
  #include "grid.h"
  #include "globalconf.h"
  #include "geom.h"
  #include "conetube.h"
  #include "pressurebc.h"
  #include "adiabaticwall.h"  
  #include "tube.h"
  #include "isentropictube.h"

%}
using std::string;
typedef std::complex<double> c;

%include "std_string.i"
%include "arma_numpy.i"
%include "std_complex.i"

%include "globalconf.h"
%include "var.h"


%include "tasystem.h"
%include "grid.h"
%include "geom.h"
%include "conetube.h"
%include "pos.h"
%include "pressurebc.h"
%include "adiabaticwall.h"  

%include "seg.h"
%include "varnr.h"
%include "tube.h"
%include "isentropictube.h"

