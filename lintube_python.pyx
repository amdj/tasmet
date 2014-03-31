import numpy as n
cimport numpy as n

cdef extern from "<vector>" namespace "std":
	cdef cppclass vector[double]:
		vector() except +
		vector(double&,double&)
		double operator[](int)
		void push_back(double)
		int size()
cdef extern from "vvar.h" namespace "vvar":
	cdef cppclass vvar:
		vvar() except +
		vector[double] getResult(unsigned)
		void setResult(vector[double],unsigned)


# cdef extern from "globalconf.h" namespace "tasystem"
# cdef extern from "tube/tube.h" namespace "tube":
# 	cdef cppclass tube:

# 		vector[double] getResult()
# 		void setResult(vector[double])
# # 		cppclass vvar:
# # 			vector[double] getResult(unsigned)
# 		vvar p
# 		vvar U
# 		vvar rho
# # 		vector[double] getUResult()
# # 		vector[double] getpResult(int)
# # 		vector[double] getUResult(int)
# 		vector[double] getError()
# 		vector[double] getx()


cdef class tube_py:
	# cdef tube* tube1
	cpdef int Nf,gp
	cpdef double freq,Up
	def __cinit__(self,int gp,int Nf,double freq,double p1):
		# self.thisl=new isentropictube(gp,Nf,1.,1.,freq,Up)
		self.Nf=Nf
		self.freq=freq
		self.Up=Up
		self.gp=gp
		self.Nf=Nf
	def getNf(self):
		return self.Nf
	def getgp(self):
		return self.gp
	def getUp(self):
		return self.Up
	def getfreq(self):
		return self.freq
	def getx(self):
		return ToArray(self.thisl.getx())
	def getResult(self):
		cdef vector[double] res=self.thisl.getResult()
		return ToArray(res)
	def getError(self):
		cdef vector[double] er=self.thisl.getError()
		return ToArray(er)
	def getpResult(self,int j):
		cdef vector[double] res=self.thisl.p.getResult(j)
		return ToArray(res)
	def getUResult(self,int j):
		cdef vector[double] res=self.thisl.U.getResult(j)
		return ToArray(res)
	def getrhoResult(self,int j):
		cdef vector[double] res=self.thisl.rho.getResult(j)
		return ToArray(res)
	def setResult(self,n.ndarray[n.float64_t,ndim=1] res):
		self.thisl.setResult(ToVec(res))
	def setpResult(self,n.ndarray[n.float64_t,ndim=1] res,int f):
		self.thisl.p.setResult(ToVec(res),f)
	def setUResult(self,n.ndarray[n.float64_t,ndim=1] res,int f):
		self.thisl.U.setResult(ToVec(res),f)
	def Solve(self):
		cdef n.ndarray[n.float64_t,ndim=1] guess=self.getResult()
		newton_raphson(self,guess,1e-3,1e-6,10,1)
	def Objfun(self,n.ndarray[n.float64_t,ndim=1] guess):
		self.setResult(guess)
		return self.getError()
	def __dealloc__(self):
		del self.thisl

# @cython.boundscheck(False)
# @cython.wraparound(False)
cdef newton_raphson(tube_py tube,n.ndarray[n.float64_t,ndim=1] guess,double reltol,double funtol,int maxiter,int verbose):
	if verbose==1:
		print "Newton-raphson iteration on function started..."
	cdef int nloop=0
	cdef int nconverg=0
	cdef double relerror,funerror
	cdef n.ndarray[n.float64_t,ndim=1] sol=guess
	cdef n.ndarray[n.float64_t,ndim=1] fx,fx2,x2,dx
	cdef int nsys=sol.shape[0]
	cdef double normx=n.linalg.norm(guess)
	cdef n.ndarray[n.float64_t,ndim=2] J
	sol=guess
	while nloop<maxiter:
		J=n.zeros((nsys,nsys))
		print "Computing zero solution"
		fx=tube.Objfun(sol)
		print "Computing jacobian..."
		for j in xrange(nsys):
			x2=sol.copy()
			unitj=n.zeros((nsys))
			unitj[j]=+0.01*normx*reltol
			x2=x2+unitj
			fx2=tube.Objfun(x2)
			dfdxj=(fx2-fx)/(x2[j]-sol[j])
			J[:,j]=dfdxj.ravel()
		if n.linalg.det(J)==0:
			print 'Determinant of Jacobian becomes zero'
		print "Solving linear system..."
		dx=-n.linalg.solve(J,fx)
		sol=sol+dx
		relerror=n.linalg.norm(dx)/n.linalg.norm(sol)
		funerror=n.linalg.norm(fx)
		nloop+=1
		nconverg+=1
		if verbose==1:
				print "Iteration: %g. Relative error: %.2e. Function error: %.2e" %(nloop,relerror,funerror)
		if nloop==maxiter-1:
			print "Warning: maximum number of iterations reached, but convergence is not reached!"
		if relerror<reltol and funerror<funtol:
			nloop=maxiter
			if verbose==1:
				print "Iteration on function converged"
				print "Relative error: %.2e. Function error: %.2e" %(relerror,funerror)
	return sol

