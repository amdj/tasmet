import numpy as n
cimport numpy as n
cdef extern from "<armadillo>" namespace "arma":
    cdef cppclass vec:
        vec(int)
        vec(double*,int,bool,bool)
        void raw_print()
        vec()
        double* memptr()
        int size()

ctypedef vec vd
ctypedef double d


cdef vd dndtovec(n.ndarray[n.float64_t,ndim=1] ndarr)
