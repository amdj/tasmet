#!/usr/bin/python

from cmath import *
import numpy as n
from copy import *
from materials import *
from var import *

class volume:
	def __init__(self,T,Nf=1):
		self.gp=[]
		self.Nf=Nf
		self.Ns=2*Nf+1
		self.T=T
		self.omg0=2*pi/T

		self.diff=n.asarray([1j*self.omg0*k for k in range(self.Nf+1)])
		self.p=var(period=T,Nf=Nf,gridpoints=1,name='Pressure',initval=101325.)
		self.T=var(period=T,Nf=Nf,gridpoints=1,name='Temperature',initval=293.15)
		self.rho=var(period=T,Nf=Nf,gridpoints=1,name='Density',initval=1.2)
		self.u=var(period=T,Nf=Nf,gridpoints=1,name='Velocity',initval=0)
		self.time=time(period=T,Nf=Nf,gridpoints=1,name='Time',initval=0)

	def __repr__(self):
		return "End volume"
	def addbc(self,bctype,val):
		lenval=len(val)
		self.bc[0:lenval]=val
	def Init(self,V,S,Tm=293.15):
		self.Tm=Tm
		self.m=air(Tm)
		self.V=V
		self.S=S
		self.bctype='pres'
		self.bc=n.zeros((self.Nf+1,),complex)
	def dft(self,var):
		Nf=self.Nf
		Ns=self.Ns
		var=n.array(var)
		ft=n.fft.rfft(var) #Fast Fourier Transform
# 		print ft
		newvar=n.zeros(Nf+1,complex)
# 		print newvar
		newvar[0:Nf+1]=ft[0:Nf+1]/(Ns)*2.
		newvar[0]*=0.5
		return newvar
	def solvec(self):
		self.sol=n.zeros(((self.Nf+1)*4,),complex)
		self.sol[0:self.Nf+1]=self.p()
		self.sol[1*(self.Nf+1):2*(self.Nf+1)]=self.rho()
		self.sol[2*(self.Nf+1):3*(self.Nf+1)]=self.T()
		self.sol[3*(self.Nf+1):4*(self.Nf+1)]=self.u()
		return self.sol
# 	def d_rho_u_dx(self,i):
# 		rho_u_ip1=

	def err(self,sol):
		sol=n.asarray(sol)
		self.p.setAdata(sol[0:self.Nf+1])
		self.u.setAdata(sol[self.Nf+1:2*(self.Nf+1)])
		self.T.setAdata(sol[2*(self.Nf+1):3*(self.Nf+1)])
		self.rho.setAdata(sol[3*(self.Nf+1):4*(self.Nf+1)])

		self.er=n.ones(((self.Nf+1)*4,),complex)

		Tt=self.T.getTdata()
		rhot=self.rho.getTdata()
		pt=self.p.getTdata()
		ut=self.u.getTdata()

		rhoe=self.dft(rhot*self.m.e_T(Tt))
		if self.bctype=='pres':
			#Eerste Nf Vgl voor pres bc
			self.er[0:self.Nf+1]=self.p()-self.bc
			#Tweede Nf Vgl voor velocity
			self.er[self.Nf+1:2*(self.Nf+1)]=self.diff*self.rho()*self.V+self.dft(ut*rhot)*self.S
			self.er[self.Nf+1]+=1j*self.u.zeroimag[0]
			#Derde Nf Vgl voor temperatuur
			self.er[2*(self.Nf+1):3*(self.Nf+1)]=self.diff*rhoe*self.V+self.dft(rhot*ut*(self.m.h_T(Tt)+0.5*ut*ut))*self.S
			self.er[2*(self.Nf+1)]=self.T()[0,0]-281.443139904
# 			self.er[2*(self.Nf+1)]=self.T()[0]-291.86

			#Vierde Nf vgl voor dichtheid
			self.er[3*(self.Nf+1):4*(self.Nf+1)]=self.rho()-self.dft(self.m.rho(pt,Tt))
		return self.er

