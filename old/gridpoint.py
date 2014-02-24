#!/uselfr/bin/python

from cmath import *
import numpy as n
from copy import *
from var import *

class gridpoint:
	def __init__(self,Nf,T,i):
# 		"""
# 		Nf: number of frequencieself other than the zero
# 		T: Time period
# 		"""
		self.i=i
		self.Nf=Nf # Num
		self.Ns=2*Nf+1
		self.period=T
		self.rho=var(T,Nf,'Density',1.2)
		self.p=var(T,Nf,'Pressure',101325.)
		self.T=var(T,Nf,'Temperature',293.15)
		self.U=var(T,Nf,'Volume flow',0.0)
		self.time=time(T,Nf,'Time',self.period)
	def product(self,varA,varB):
		return varA.dft(varA.getTdata()*varB.getTdata())
	def __repr__(self):
		return "Gridpoint at position %g" %self.i

	def setT(self,T):
		self.T=T

