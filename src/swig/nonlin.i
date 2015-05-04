%module(docstring="TATwente nonlinear thermoacoustic code") TATwente
%{
  #define PY_ARRAY_UNIQUE_SYMBOL npy_array
  #define SWIG_FILE_WITH_INIT

  // Conversion between numpy and Armadillo
  #include "arma_numpy.h"

  #include "settracer.h"
  #include "constants.h"

  // My exceptions
  #include "exception.h"

  // Global config
  #include "globalconf.h"
  #include "var.h"

  // Build a system
  #include "tasystem.h"
  #include "enginesystem.h"

  #include "boundarylayer.h"
  #include "geom.h"
  #include "conetube.h"

  // For Segments and connectors
  #include "segconbase.h"

  // Connectors
  #include "connector.h"
  #include "tubebc.h"
  #include "tubeconnector.h"

  #include "pressurebc.h"
  #include "adiabaticwall.h"  
  #include "isotwall.h"  
  // Segments
  #include "seg.h"

  // Tubes
  #include "grid.h"
  #include "tube.h"
  #include "isentropictube.h"
  #include "laminarduct.h"
  #include "hopkinslaminarduct.h"

  // Solver
  #include "solverconfiguration.h"
  #include "solver.h"

  // A small wrapper for this function, as TRACERNAME is not
  // substituted by its macro value in SWIG.
  inline void setTATwenteTracer(int t) {
    tracer::setTracer<TRACERNAME>(t);
  }

  %}
using std::string;
typedef std::complex<double> c;

// To change the tracing noisyness
void setTATwenteTracer(int);


%include "std_string.i"
%include "std_complex.i"

 // Conversion between numpy and Armadillo
%include "arma_numpy.i"

 // My compile-time constants

%include "constants.h"

// My exceptions
%include "my_exceptions.i"

 // Global config
%feature("autodoc","3");
%include "globalconf.h"
%include "var.h"


 // For Segments and connectors
%include "segconbase.h"
%include "phaseconstraint.h"
 // Connectors
%include "connector.h"
%include "tubebc.h"
%include "tubeconnector.h"

// Other connectors
%include "pressurebc.h"
%include "adiabaticwall.h"  
%include "isotwall.h"  
 // Segments
%include "seg.h"

 // Tubes
%include "grid.h"
%include "boundarylayer.h"
%include "geom.h"
 // Geom instances
%include "conetube.h"
%include "vertplates.h"
 // Tubes
%include "tube.h"
%include "isentropictube.h"
%include "laminarduct.h"
%include "hopkinslaminarduct.h"

 // System
%include "tasystem.h"
%include "enginesystem.h"
 // Solver
%include "solverconfiguration.h"
%include "solver.h"
