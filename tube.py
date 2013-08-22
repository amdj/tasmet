from cmath import *
from math import *
import numpy as n
from var import *

from materials import air
from testbc import *

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
		newvar[0,Nf + 1:Ns] = (var[0,Nf + 1:0:-1]*Ns*0.5).conjugate()
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

class tube(segment,testing,testbc):
	def __init__(self, period, Nf, L, S, gridpoints=10,x=None):
		self.segInit()
		self.gridpoints=gridpoints
		self.Nf = Nf
		self.Ns = 2 * Nf + 1
		self.period = period
		self.omg0=2*pi/self.period
		self.testInit()
		self.omega=n.zeros((1,Nf+1),float)
		for i in xrange(Nf+1):
			self.omega[0,i]=self.omg0*i
		self.m=air()
		self.rho=var(period,Nf,gridpoints,name="Density",initval=1.2)
		self.U=var(period,Nf,gridpoints,name="Volume flow",initval=1e-6)
		self.p=var(period,Nf,gridpoints,name="Pressure",initval=100000)
		self.T=var(period,Nf,gridpoints,name="Temperature",initval=300)

		self.mf=var(period,Nf,gridpoints,name="Mass flow",initval=1)
		self.momf=var(period,Nf,gridpoints,name="Momentum flow",initval=1)
		self.Etot=var(period,Nf,gridpoints,name="Total energy",initval=self.m.e_T(293.15))
		self.Hf=var(period,Nf,gridpoints,name="Total enthalpy flow",initval=0)

		if x is None:
			self.x=n.linspace(0,L,gridpoints)
			self.Sf=S*n.ones((self.gridpoints,self.Nf+1),float)
			self.SfT=S*n.ones((self.gridpoints,self.Ns),float)
		self.ndofs=gridpoints*4*self.Ns
		self.iter=0
		self.Update()
		#for i in xrange(gridpoints):
		#	Si=S
		#	self.gp.append(gridpoint(self.Nf,T,i))
	def __repr__(self):
		return "Tube class"
	def blocks(self,i):
		blk=self.Nf+1
		rhoblk=[4*i*blk,4*i*blk+blk]
		Ublk=[4*i*blk+blk,4*i*blk+2*blk]
		pblk=[4*i*blk+2*blk,4*i*blk+3*blk]
		Tblk=[4*i*blk+3*blk,4*i*blk+4*blk]
		return rhoblk,Ublk,pblk,Tblk
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
		newT=n.zeros(self.rho.shape,complex)
		#This CAN BE PARALELLIZED
		for i in xrange(self.gridpoints):
			rhoblk,Ublk,pblk,Tblk=self.blocks(i)
			newrho[i]=cres[rhoblk[0]:rhoblk[1]]
			newU[i]=cres[Ublk[0]:Ublk[1]]
			newp[i]=cres[pblk[0]:pblk[1]]
			newT[i]=cres[Tblk[0]:Tblk[1]]
		self.rho.setAdata(newrho)
		self.U.setAdata(newU)
		self.T.setAdata(newT)
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
		T=self.T.getAdata()
		cres=n.zeros((4*(self.Nf+1)*self.gridpoints,),complex)
		#cdef unsigned int i
		for i in xrange(self.gridpoints):
			rhoblk,Ublk,pblk,Tblk=self.blocks(i)
			cres[rhoblk[0]:rhoblk[1]]=rho[i].ravel()
			cres[Ublk[0]:Ublk[1]]=U[i].ravel()
			cres[pblk[0]:pblk[1]]=p[i].ravel()
			cres[Tblk[0]:Tblk[1]]=T[i].ravel()
# 			#=rho[0].ravel()
		res=self.complextofull(cres)
		return res

	def Update(self):
		u=self.U.getTdata()/self.SfT
		self.mf.setTdata(self.rho.getTdata()*self.U.getTdata())
		self.momf.setTdata(self.mf.getTdata()*u)
		# FIX THIS!!!
		Ekin=0.5*u*u
		self.Etot.setTdata(self.rho.getTdata()*(self.m.e_T(self.T.getTdata())+Ekin))
		self.Hf.setTdata(self.mf.getTdata()*(self.m.h_T(self.T.getTdata())+Ekin))
	def setImagzero(self,error):
		for i in xrange(0,self.gridpoints):
			rhoblk,Ublk,pblk,Tblk=self.blocks(i)
			error[rhoblk[0]]+=1j*self.rho.zeroimag[i]
			error[Tblk[0]]+=1j*self.T.zeroimag[i]
			#Set imaginary part of mean pressure to zero
			error[pblk[0]]+=1j*self.p.zeroimag[i]
			error[Ublk[0]]+=1j*self.U.zeroimag[i]
		return error
	def innerError(self,error):
		#Can be parallellized:3
		for i in xrange(1,self.gridpoints-1):
			rhoblk,Ublk,pblk,Tblk=self.blocks(i)
			#Implementation of continuity equation
			Vi=self.Sf[i]*0.5*(self.x[i+1]-self.x[i-1])
			dmdt=1j*self.omega*self.rho()[i]*Vi

			rhoudotn=0.5*(self.mf()[i+1]-self.mf()[i-1])
			error[rhoblk[0]:rhoblk[1]]=dmdt+rhoudotn

			#Implementation of momentum equation
			dmfdt=1j*self.omega*self.mf()[i]*Vi
			pdotn=0.5*(self.Sf[i+1]*self.p()[i+1]-self.Sf[i-1]*self.p()[i-1])

			rhouudotn=0.5*(self.momf()[i+1]-self.momf()[i-1])
			#FIXTHIS
			error[Ublk[0]:Ublk[1]]=dmfdt+pdotn+rhouudotn

			#Implementation of energy equation
			dEtotdt=1j*self.omega*self.Etot()[i]*Vi
			error[Tblk[0]:Tblk[1]]=dEtotdt+0.5*(self.Hf()[i+1]-self.Hf()[i-1])
			error[Tblk[0]]=self.T()[i,0]-293.15

			#Implementation of equation of state
			#print self.dft(self.m.p_rho_T(self.rho.getTdata()[i], self.T.getTdata()[i]))
			error[pblk[0]:pblk[1]]=self.p()[i]-self.dft(self.m.p_rho_T(self.rho.getTdata()[i], self.T.getTdata()[i]))

		return error

	def getError(self):
		error=n.zeros((4*(self.Nf+1)*self.gridpoints,),complex)
		error=self.innerError(error)
		error=self.leftBc3(error)
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
