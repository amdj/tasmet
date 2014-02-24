#!/usr/bin/python
import numpy as n
from numpy import exp,pi
from lintube_python import *
import cPickle as pickle
from tmtubes.linear import *
class storeobj:
	def __init__(self,t):
		self.U=[]
		self.p=[]
		self.rho=[]
		self.Nf=t.getNf()
		self.gp=t.getgp()
		self.Up=t.getUp()
		self.freq=t.getfreq()
		self.p.append(t.getpResult(0))
		self.U.append(t.getUResult(0))
		self.rho.append(t.getrhoResult(0))	
		for i in range(1,self.Nf+1):
			self.U.append(t.getUResult(2*i-1)+1j*t.getUResult(2*i))
			self.rho.append(t.getrhoResult(2*i-1)+1j*t.getrhoResult(2*i))
			self.p.append(t.getpResult(2*i-1)+1j*t.getpResult(2*i))
		self.x=t.getx()
		self.filename='gp%g-Nf%g-Up%g-freq%g' %(self.gp,self.Nf,self.Up,self.freq)
		ofile=open(self.filename+'.dat','w')
		pickle.dump(self,ofile)
		ofile.close()
		t1=ISOTDUCT()
		t1.setgrid(L=1,gridpoints=self.gp)
		t1.setnodes(0,1)
		t1.creategeom(cshape=None,S=1.)
		t1.Init(294.20731707317077)
		s1=TAsystem([t1])
		s1.acousticSystem.addbc(0,'volu',self.Up)
		s1.acousticSystem.addbc(1,'volu',0)
		s1.Init('air',freq=self.freq)
		s1.Acfreqresponse(self.freq)
		self.lt=s1.segs[0]
		self.linx=self.lt.x
	def linUt(self,t):
		return (self.lt.U1*exp(1j*2*pi*self.freq*t)).real
	def linpt(self,t):
		return (self.lt.p1*exp(1j*2*pi*self.freq*t)).real
	def Ut(self,t):
		res=n.zeros(self.x.shape[0],)
		for i in range(self.Nf):
			res=res+(self.U[i]*exp(1j*2*pi*self.freq*i*t)).real
		return res
	def rhot(self,t):
		res=n.zeros(self.x.shape[0],)
		for i in range(self.Nf):
			res=res+(self.rho[i]*exp(1j*2*pi*self.freq*i*t)).real
		return res
	def pt(self,t):
		res=n.zeros(self.x.shape[0],)
		for i in range(self.Nf):
			res=res+(self.p[i]*exp(1j*2*pi*self.freq*i*t)).real
		return res
