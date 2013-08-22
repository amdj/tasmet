from cmath import *
from math import *
import numpy as n
from var import *

from materials import air
from testlinbc import *

class segment:
	number=0
	def segInit(self):
		self.number=segment.number
		segment.number+=1
	def idft1long(self,  var):
		Nf = self.Nf
		Ns = self.Ns
		newvar = n.zeros((1,self.Ns), complex)
		newvar[0,1:Nf + 1] = var[0,1:Nf + 1] * Ns * 0.5
		newvar[0,Nf + 1:Ns] = (var[0,Nf + 1:0:-1] * Ns * 0.5).conjugate()
		#print newvar
		#print n.fft.ifft(newvar)
		return n.fft.ifft(newvar).real
	def dft(self, var):
		Nf = self.Nf
		Ns = self.Ns
		if len(var.shape)==1:
			rows=1
			newvar = n.zeros((rows,self.Nf+1), complex)
			ft = n.fft.rfft(var)  # Fast Fourier Transform
			newvar[0:Nf + 1] = ft[0:Nf + 1]/(Ns)*2.
			newvar[0]*= 0.5
		else:
			rows=var.shape[0]
			newvar = n.zeros((rows,self.Nf+1), complex)
			ft = n.fft.rfft(var)  # Fast Fourier Transform
			newvar[:,0:Nf + 1] = ft[:,0:Nf + 1]/(Ns)*2.
			newvar[:,0]*= 0.5
		return newvar

	def complextofull(self,complexval):
		length=complexval.shape[0]
		full=n.zeros((2*length,),float)
		full[0:2*length:2]=complexval.real[:]
		full[1:2*length+1:2]=complexval.imag[:]
		return full
	def fulltocomplex(self,full):
		length=full.shape[0]/2
		#print "full shape:",full.shape
		#cdef n.ndarray[dtype=n.complex_t,ndim=1] complexval=n.zeros((length,),complex)
		complexval=n.zeros((length,),complex)
		complexval[0:length]=full[0:2*length:2]
		complexval[0:length]=complexval[0:length]+1j*full[1:2*length+1:2]
		return complexval

class lintube(segment,testlinbc):
	def __init__(self, period, Nf, L, S, gridpoints=10,x=None):
		self.segInit()
		self.gridpoints=gridpoints
		self.Nf = Nf
		self.Ns = 2 * Nf + 1
		self.period = period
		self.omg0=2*pi/self.period
		self.omega=n.zeros((1,Nf+1),float)
		for i in xrange(Nf+1):
			self.omega[0,i]=self.omg0*i
		self.m=air()
		self.rho=var(period,Nf,gridpoints,name="Density",initval=1.2)
		self.U=var(period,Nf,gridpoints,name="Volume flow",initval=1e-6)
		self.p=var(period,Nf,gridpoints,name="Pressure",initval=100000)

		if x is None:
			self.x=n.linspace(0,L,gridpoints)
			self.Sf=S*n.ones((self.gridpoints,self.Nf+1),float)
			self.SfT=S*n.ones((self.gridpoints,self.Ns),float)
		self.ndofs=gridpoints*3*self.Ns
		self.iter=0
		self.Update()
		#for i in xrange(gridpoints):
		#	Si=S
		#	self.gp.append(gridpoint(self.Nf,T,i))
	def __repr__(self):
		return "Tube class"
	def blocks(self,i):
		blk=self.Nf+1
		rhoblk=[3*i*blk,3*i*blk+blk]
		Ublk=[3*i*blk+blk,3*i*blk+2*blk]
		pblk=[3*i*blk+2*blk,3*i*blk+3*blk]
		return rhoblk,Ublk,pblk
	def setResult(self, result):
		#Every block is a grid point
		#In every block, first the density dofs
		#Then Volume flow dofs
		#Then pressure dofs
		#Then Temperature dofs

		#First we convert from double length reals to complex-valued
		cres=self.fulltocomplex(result)
# 		print cres
		newrho=n.zeros(self.rho.shape,complex)
		newU=n.zeros(self.rho.shape,complex)
		newp=n.zeros(self.rho.shape,complex)

		#This CAN BE PARALELLIZED
		for i in xrange(self.gridpoints):
			rhoblk,Ublk,pblk=self.blocks(i)
			newrho[i]=cres[rhoblk[0]:rhoblk[1]]
			newU[i]=cres[Ublk[0]:Ublk[1]]
			newp[i]=cres[pblk[0]:pblk[1]]
		self.rho.setAdata(newrho)
		self.U.setAdata(newU)
		self.p.setAdata(newp)
		#And update others...
		self.Update()
	def getResult(self):
		#In every block, first the density dofs
		#Then Volume flow dofs
		#Then pressure dofs
		#Then Temperature dofs
		rho=self.rho.getAdata()
		U=self.U.getAdata()
		p=self.p.getAdata()
		cres=n.zeros((3*(self.Nf+1)*self.gridpoints,),complex)
		#cdef unsigned int i
		for i in xrange(self.gridpoints):
			rhoblk,Ublk,pblk=self.blocks(i)
			cres[rhoblk[0]:rhoblk[1]]=rho[i].ravel()
			cres[Ublk[0]:Ublk[1]]=U[i].ravel()
			cres[pblk[0]:pblk[1]]=p[i].ravel()
# 			#=rho[0].ravel()
		res=self.complextofull(cres)
		return res
	def Update(self):
		pass
	def setImagzero(self,error):
		for i in xrange(0,self.gridpoints):
			rhoblk,Ublk,pblk=self.blocks(i)
			error[rhoblk[0]]+=1j*self.rho.zeroimag[i]
			#Set imaginary part of mean pressure to zero
			error[pblk[0]]+=1j*self.p.zeroimag[i]
			error[Ublk[0]]+=1j*self.U.zeroimag[i]
		return error
	def innerError(self,error):
		#Can be parallellized:3
		for i in xrange(1,self.gridpoints-1):
			rhoblk,Ublk,pblk=self.blocks(i)
			dUdt=1j*self.omega*self.U()[i]
			drhodt=1j*self.omega*self.rho()[i]
			Udotn=(self.U()[i+1]-self.U()[i-1])
			Vi=self.x[i+1]-self.x[i-1]
			error[rhoblk[0]:rhoblk[1]]=Vi*drhodt+1.2*Udotn
			error[rhoblk[0]]=self.rho()[i,0]-1.2
			#Implementation of momentum equation
			pdotn=(self.p()[i+1]-self.p()[i-1])
			error[Ublk[0]:Ublk[1]]=1.2*Vi*dUdt+pdotn
			error[Ublk[0]]=self.U()[i,0]
			error[pblk[0]:pblk[1]]=self.p()[i]-287*1.4*293.15*self.rho()[i]
			error[pblk[0]]=self.p()[i,0]-101325.
		return error
	def getError(self):
		error=n.zeros((3*(self.Nf+1)*self.gridpoints,),complex)
		error=self.innerError(error)
		error=self.leftBc(error)
		# Boundary conditions, manually
		error=self.rightBc(error)
		error=self.setImagzero(error)
		rerror=self.complextofull(error)
		self._error=rerror
		return rerror
	def Solve(self,result):
		#print "Iteration:" ,self.iter
		self.iter+=1
		#print "Result shape from solver:" , result.shape
		self.setResult(result)
		error=self.getError()
		#print "Error shape:" ,error.shape
		#print self.fulltocomplex(error)
		return error
#
# 		#Using adiabatic gas law
# 		T0=293.15
# 		p0=101325.
#
# 		#intcp_overTdT=self.dft(self.m.intcp_over_TdT(T0, self.T.getTdata()[i]))
# 		#int1_overpdp=self.dft(n.log(self.p.getTdata()[i]/p0))
#
# 		#error[Tblk[0]:Tblk[1],0]=intcp_overTdT-int1_overpdp/self.m.Rs
# 		#error[Tblk[0],0]+=1j*self.T.zeroimag[i]
