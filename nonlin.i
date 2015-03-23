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

  #include "solver.h"
%}
using std::string;
typedef std::complex<double> c;

%include "std_string.i"
%include "std_complex.i"
%include "arma_numpy.i"

// My exceptions
%include "my_exceptions.i"

 // Global config
%include "globalconf.h"
%include "var.h"


%include "tasystem.h"
%include "grid.h"
%include "geom.h"
%include "conetube.h"
%include "pos.h"

 // Connectors
%include "connector.h"
%include "tubebc.h"

%include "pressurebc.h"
%include "adiabaticwall.h"  

 // Segments
%include "seg.h"
%include "varnr.h"
%include "tube.h"
%include "isentropictube.h"

 // Solver
%include "solver.h"
