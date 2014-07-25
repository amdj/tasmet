include "pytube.pxi"

cdef extern from "fubini.h" namespace "":
    Solver* Fubini(us gp,us Nf,d freq,d L,d S,vd p1,int loglevel,d kappa)
cdef extern from "threetubes.h" namespace "":
    Solver* ThreeTubes(us gp,us Nf,d freq,d L,d S,vd p1,int loglevel,d kappa)


#Pytube compatible with paper_1 Fubini code
    
cdef class pytube:
    cdef Solver* sol
    cdef Tube* tube[5]
    def __cinit__(self,us gp,us Nf,d freq,d L,d S,d T0,d p0,n.ndarray[n.float64_t,ndim=1] p1, cshape,int loglevel,d kappa,case="fubini"):
        if case is "fubini":
            print('Case Fubini')
            self.sol=Fubini(gp,Nf,freq,L,S,dndtovec(p1),loglevel,kappa)
            self.tube[0]=<Tube*> self.sol[0].sys[0].getSeg(0)
        elif case is "threetubes":
            self.sol=ThreeTubes(gp,Nf,freq,L,S,dndtovec(p1),loglevel,kappa)
            self.tube[0]=<Tube*> self.sol[0].sys[0].getSeg(0)
            self.tube[1]=<Tube*> self.sol[0].sys[0].getSeg(1)
        else:
            print('Warning: no valid case selected! Tried was: %s' %case)

    def show(self):
        self.sol[0].sys[0].show()    
    def __dealloc__(self):
        del self.sol
    cpdef getNf(self):
        return self.Nf
    cpdef getgp(self):
        return self.gp
    cpdef getx(self,i=0):
        return dvectond(self.tube[i].geom.vx)
    cpdef Error(self):
        return dvectond(self.sol[0].sys[0].Error())
    cpdef GetRes(self):
        return dvectond(self.sol[0].sys[0].GetRes())
    cpdef DoIter(self,d relaxfac):
        self.sol[0].DoIter(relaxfac)
    cpdef SetRes(self,n.ndarray[n.float64_t,ndim=1] res):
        self.sol[0].sys[0].SetRes(dndtovec(res))
    cpdef GetResVar(self,_type,freqnr,i=0):
        if _type=='pres':
            return dvectond(self.tube[i].GetResAt(3,freqnr))
        elif _type=='rho':
            return dvectond(self.tube[i].GetResAt(0,freqnr))
        elif _type=='temp':
            return dvectond(self.tube[i].GetResAt(2,freqnr))
        elif _type=='volu':
            return dvectond(self.tube[i].GetResAt(1,freqnr))
        else:
            return None
    # def getResult(self):
    #     cdef vector[double] res=self.thisl.getRes()
    #     return ToArray(res)
    # def getError(self):
    #     cdef vector[double] er=self.thisl.getErr()
    #     return ToArray(er)
    # def getpResult(self,int j):
    #     cdef vector[double] res=self.thisl.p.getRes(j)
    #     return ToArray(res)
    # def getUResult(self,int j):
    #     cdef vector[double] res=self.thisl.U.getRes(j)
    #     return ToArray(res)
    # def getrhoResult(self,int j):
    #     cdef vector[double] res=self.thisl.rho.getRes(j)
    #     return ToArray(res)
    # def setResult(self,n.ndarray[n.float64_t,ndim=1] res):
    #     self.thisl.setResult(ToVec(res))
    # def setpResult(self,n.ndarray[n.float64_t,ndim=1] res,int f):
    #     self.thisl.p.setResult(ToVec(res),f)
    # def setUResult(self,n.ndarray[n.float64_t,ndim=1] res,int f):
    #     self.thisl.U.setResult(ToVec(res),f)
    # def Objfun(self,n.ndarray[n.float64_t,ndim=1] guess):
    #     self.setResult(guess)
    #     return self.getError()
        # del self.tube1

