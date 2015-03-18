%module nonlin
%{
  #define PY_ARRAY_UNIQUE_SYMBOL npy_array
  #include "globalconf.h"
%}
%include "std_string.i"
%include "arma_numpy.i"
%include "std_complex.i"
typedef std::string string;
typedef std::complex<double> c;

