%module(docstring="TASMET nonlinear thermoacoustic code") tasmet
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
  #include "connectorvolume.h"
  // Connectors
  #include "connector.h"
  #include "ductbc.h"
  #include "ductconnector.h"
  #include "ductpistonconnector.h"

  #include "pressurebc.h"
  #include "impedancebc.h"
  #include "velocitybc.h"
  #include "adiabaticwall.h"  
  #include "isotwall.h"  

  // Segments
  #include "seg.h"
  #include "piston.h"
  #include "mechbc.h"

  // Ducts
  #include "grid.h"
  #include "duct.h"
  #include "cell.h"
  #include "isentropictube.h"
  #include "laminarduct.h"
  #include "stack.h"

  // Solver
  #include "solverconfiguration.h"
  #include "errorvals.h"
  #include "solprogress.h"
  #include "solver.h"

  // A small wrapper for this function, as TRACERNAME is not
  // substituted by its macro value in SWIG.
  inline void setTASMETTracer(int t) {
    tracer::setTracer<TRACERNAME>(t);
  }

  %}

// Grab a Python function object as a Python object.
%typemap(in) PyObject* pyfunc {
  if (!PyCallable_Check($input)) {
      PyErr_SetString(PyExc_TypeError, "Need a callable object!");
      return NULL;
  }
  $1 = $input;
}

using std::string;
typedef std::complex<double> c;

void clearConsole();
// To change the tracing noisyness
void setTASMETTracer(int);


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
%include "ductbc.h"
%include "ductconnector.h"
%include "ductpistonconnector.h"

// Other connectors
%include "pressurebc.h"
%include "adiabaticwall.h"  
%include "isotwall.h"  
%include "impedancebc.h"
%include "velocitybc.h"
%include "mechbc.h"
 // Segments
%include "seg.h"
%include "piston.h"
%include "connectorvolume.h"
 // Ducts
%include "grid.h"
%include "boundarylayer.h"
%include "geom.h"
 // Geom instances
%include "conetube.h"
%include "vertplates.h"
 // Ducts
%include "cell.h"
%include "duct.h"
%include "isentropictube.h"
%include "laminarduct.h"
%include "ductwithsolid.h"
%include "stack.h"

 // System
%include "tasystem.h"
%include "enginesystem.h"
 // Solver
%include "errorvals.h"
%include "solprogress.h"
%include "solverconfiguration.h"
%include "solver.h"

%pythoncode{
import sys,os

sys.path.append(os.path.join(os.path.dirname(__file__),'src/gui'))
import tasmet_main

if __name__=='__main__':
    tasmet_main.run()
  
}
