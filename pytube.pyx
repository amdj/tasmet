include "pytube.pxi"

cdef extern from "fubini.h" namespace "":
    Solver* Fubini(us gp,us Nf,d freq,d L,d S,c p1,int loglevel,d kappa)

cdef class pytube:
    cdef Geom* geom1
    cdef Globalconf* gc
    cdef Solver* sol

    def __cinit__(self,us gp,us Nf,d freq,d L,d S,d T0,d p0,n.ndarray[n.float64_t,ndim=1] p1,string cshape,int loglevel,d kappa):
        # self.thisl=new isentropictube(gp,Nf,1.,1.,freq,Up)
        # print "New tube initialized"
        initlog(loglevel)
        self.Nf=Nf
        self.freq=freq
        # print "freq:",freq
        self.gp=gp
        self.Nf=Nf
        # print "p0:",p0
        cdef d phi=1.0
        cdef d R=n.sqrt(S/n.pi)
        cdef d rh=S/(2*n.pi*R)
        cdef d Mach=abs(p1[1]**2+p1[2]**2)/p0
        # print "rh:",rh
        print "cshape:",cshape
        cdef d Mass=0
        cdef d dx=L/float(gp)
        # print "dx:", dx
        # print "kappa:",kappa
        # cdef d S
        # print p1
        cdef Geom* geom1=new Geom(gp,L,S,phi,rh,cshape)
        self.tube1=new Tube(self.geom1[0])
        del self.geom1

        self.gc=new Globalconf(Nf,freq,"air",T0,p0,Mach,S,dx,Mass,kappa)
        self.pL=new var(self.gc[0])
        if(Nf>0):
            self.pL.set(dndtovec(p1))
        self.lp=new LeftPressure(self.tube1[0],self.pL[0])
        del self.pL
        self.tube1.setLeftbc(self.lp)     #tube1 owns the bc vertex, so do not delete the memory!


    cpdef setRightImpedance(self,n.ndarray[n.float64_t,ndim=1] Z):
        print "setRightImpedance called."
        cdef vd Zvec=dndtovec(Z)
        self.ri=new RightImpedance(self.tube1.get()[0],Zvec)
        self.tube1.setRightbc(self.ri)
        self.sys=new TAsystem(self.gc[0])
        self.sys.addseg(<Seg&> self.tube1)
        self.sol=new Solver(self.sys[0])
    def __dealloc__(self):
        del self.sol
        del self.sys
        del self.tube1
        del self.gc
        del self.geom1
    cpdef getNf(self):
        return self.Nf
    cpdef getgp(self):
        return self.gp
    cpdef getx(self):
        return dvectond(self.tube1[0].geom.vx)
    cpdef Error(self):
        return dvectond(*self.sys.Error())
    cpdef GetRes(self):
        return dvectond(*self.sys.GetRes())
    cpdef DoIter(self,d relaxfac):
        self.sol.DoIter(relaxfac)
    cpdef SetRes(self,n.ndarray[n.float64_t,ndim=1] res):
        self.tube1.SetRes(dndtovec(res))
    cpdef GetResVar(self,_type,freqnr):
        if _type=='pres':
            return dvectond(*self.tube1.GetResAt(3,freqnr))
        elif _type=='rho':
            return dvectond(*self.tube1.GetResAt(0,freqnr))
        elif _type=='temp':
            return dvectond(*self.tube1.GetResAt(2,freqnr))
        elif _type=='volu':
            return dvectond(*self.tube1.GetResAt(1,freqnr))
        
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

# @cython.boundscheck(False)
# @cython.wraparound(False)
# cdef newton_raphson(tube_py tube,n.ndarray[n.float64_t,ndim=1] guess,double reltol,double funtol,int maxiter,int verbose):
#   if verbose==1:
#       print "Newton-raphson iteration on function started..."
#   cdef int nloop=0
#   cdef int nconverg=0
#   cdef double relerror,funerror
#   cdef n.ndarray[n.float64_t,ndim=1] sol=guess
#   cdef n.ndarray[n.float64_t,ndim=1] fx,fx2,x2,dx
#   cdef int nsys=sol.shape[0]
#   cdef double normx=n.linalg.norm(guess)
#   cdef n.ndarray[n.float64_t,ndim=2] J
#   sol=guess
#   while nloop<maxiter:
#       J=n.zeros((nsys,nsys))
#       print "Computing zero solution"
#       fx=tube.Objfun(sol)
#       print "Computing jacobian..."
#       for j in xrange(nsys):
#           x2=sol.copy()
#           unitj=n.zeros((nsys))
#           unitj[j]=+0.01*normx*reltol
#           x2=x2+unitj
#           fx2=tube.Objfun(x2)
#           dfdxj=(fx2-fx)/(x2[j]-sol[j])
#           J[:,j]=dfdxj.ravel()
#       if n.linalg.det(J)==0:
#           print 'Determinant of Jacobian becomes zero'
#       print "Solving linear system..."
#       dx=-n.linalg.solve(J,fx)
#       sol=sol+dx
#       relerror=n.linalg.norm(dx)/n.linalg.norm(sol)
#       funerror=n.linalg.norm(fx)
#       nloop+=1
#       nconverg+=1
#       if verbose==1:
#               print "Iteration: %g. Relative error: %.2e. Function error: %.2e" %(nloop,relerror,funerror)
#       if nloop==maxiter-1:
#           print "Warning: maximum number of iterations reached, but convergence is not reached!"
#       if relerror<reltol and funerror<funtol:
#           nloop=maxiter
#           if verbose==1:
#               print "Iteration on function converged"
#               print "Relative error: %.2e. Function error: %.2e" %(relerror,funerror)
#   return sol
