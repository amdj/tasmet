#!/usr/bin/python

from cmath import *
import numpy as n
from copy import *

class gridpoint:
	def __init__(s,Nf,T):
# 		"""
# 		Nf: number of frequencies other than the zero
# 		T: Time period
# 		"""
		s.T=float(T)
		s.Nf=Nf # Num
		s.Ns=2*Nf+1

		s.omg0=2*pi/T
		s.rho=var(s,'Density',1.2)
# 		s.pres=var(s,'Pressure',101325.)
# 		s.temp=var(s,'Temperature',293.15)
		s.u=var(s,'Velocity',0.0)
		s.time=time(s,'Time',s.T)
	def product(s,varA,varB):
		return varA.dft(varA.getTdata()*varB.getTdata())

	def setT(s,T):
		s.T=T

class var:
	def __init__(s,parent,name,initval):
		s.name=name
		s.T=parent.T
		s.Nf=parent.Nf # Number of frequencies
		s.Ns=2*s.Nf+1
		s.omg0=2*pi/s.T

		s.tdata=initval*n.ones((s.Ns,))
		s.adata=s.dft(copy(s.tdata))
	def __call__(s):
		return s.adata()
	def setAdata(s,var):
		lenvar=len(var)
		print s.adata.shape
		print "blabla"
		s.adata[0:lenvar]=n.asarray(var)
		s.tdata=s.idft(deepcopy(n.asarray(var)))
	def setTdata(s,var):
		s.tdata=n.asarray(var)
		s.adata=s.dft(deepcopy(s.tdata))
	
	def getTdata(s):
		print "blabla"
		return s.tdata
	def getAdata(s):
		print "blabla"
		return s.adata

	def show(s):
		print "Name:",s.name
		print "Amplitudes:",s.listAdata()

	def listTdata(s):
		returnval=s.tdata.ravel().tolist()
		return returnval

	def listAdata(s):
		return s.adata.ravel().tolist()
	
	def getProfile(s,periods=1,points=None):
		if points==None:
			points=s.Ns
		times=n.arange(0,s.T*periods+s.T*periods/points,s.T*periods/points)
		varlist=n.zeros(len(times.tolist()),complex)
		alist=s.listAdata()
		for freqnr in range(s.Nf+1):
			varlist+=alist[freqnr]*n.exp(1j*s.omg0*freqnr*times)
		return varlist.real.tolist()

	def dft(s,var):
		var=n.array(var)
		ft=2*n.fft.rfft(var)/s.Ns #Fast Fourier Transform
		ft[0]*=0.5 #Multiply steady value by half
		return ft	

	def idft(s,var):
		var=n.asarray(var)
# 		var[0]*=2.
		var*=s.Ns/2.
# 		s.Nf=s.s.Nff #Shorten variable
		print var.shape[0]
		newvar= n.asarray(n.zeros((var.shape[0]*2-1),complex ))
		newvar[s.Nf-1:]=var[:s.Nf+1]
		newvar[0:s.Nf]=var[s.Nf:0:-1].conjugate()
		print newvar
		return n.fft.ifft(newvar).real
class time(var):
	def __init__(s,parent,name,initval):
		s.name=name
		s.T=parent.T
		s.Nf=parent.Nf # Number of frequencies
		s.Ns=2*s.Nf+1
		s.omg0=2*pi/s.T

		s.tdata=n.arange(0,s.T,s.T/s.Ns)
		s.adata=s.dft(copy(s.tdata))

	def getProfile(s,periods=1,points=None):
		if points==None:
			points=s.Ns
		times=n.arange(0,s.T*periods+s.T*periods/points,s.T*periods/points)
		return times
