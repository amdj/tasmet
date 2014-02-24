'''
Created on Jul 30, 2013

@author: anne
'''
from math import *
import numpy as n



class var:
	def __init__(self, period, Nf, gridpoints, name, initval):
		self.name = name
		self.T = period
		self.Nf = Nf  # Number of frequencies
		self.Ns = 2 * self.Nf + 1
		self.omg0 = 2 * pi / self.T
		self.gridpoints = gridpoints
		self.tdata = initval * n.ones((gridpoints,self.Ns),float)
		self.adata = self.dft(n.copy(self.tdata))
		self.zeroimag=n.zeros((gridpoints,),float)
		self.shape=(gridpoints,self.Nf+1)
		self.Tshape=(gridpoints,2*self.Nf+1)

	def __call__(self):
		return self.adata
	def setAdata(self, var):
		lenvar = len(var)
		self.adata[:,:] = n.asarray(var)
		self.tdata = self.idft(n.copy(self.adata))
	def setTdata(self, var):
		self.tdata = var
		self.adata = self.dft(n.copy(self.tdata))
	def getTdata(self):
		return self.tdata
	def getAdata(self):
		return self.adata

	def show(self):
		print "Name:", self.name
		print "Amplitudes:", self.listAdata()

	def listTdata(self):
		returnval = self.tdata.ravel().tolist()
		return returnval

	def listAdata(self):
		return self.adata.ravel().tolist()

	def getProfile(self, points=None, periods=1):
		if points == None:
				points = self.Ns
		times = n.arange(0, self.T * periods + self.T * periods / points, self.T * periods / points)
		varlist = n.zeros(len(times.tolist()), complex)
		alist = self.listAdata()
		for freqnr in range(self.Nf + 1):
				varlist += alist[freqnr] * n.exp(1j * self.omg0 * freqnr * times)
		return varlist.real.tolist()

	def dft(self, var):
		Nf = self.Nf
		Ns = self.Ns
		ft = n.fft.rfft(var)  # Fast Fourier Transform
		newvar = n.zeros((self.gridpoints,self.Nf+1), complex)
		newvar[:,0:Nf + 1] = ft[:,0:Nf + 1]/(Ns)*2.
		newvar[:,0]*= 0.5
		return newvar

	def idft(self, var):
		Nf = self.Nf
		Ns = self.Ns
		newvar = n.asarray(n.zeros((self.gridpoints,self.Ns), complex))
		newvar[:,0] = var[:,0] * Ns
		newvar[:,1:Nf + 1] = var[:,1:Nf + 1] * Ns * 0.5
		newvar[:,Nf + 1:Ns] = (var[:,Nf + 1:0:-1] * Ns * 0.5).conjugate()
		#print newvar
		#print n.fft.ifft(newvar)
		ifft=n.fft.ifft(newvar)
		self.zeroimag=ifft[:,0].imag
		return ifft.real
class time(var):
	def getProfile(self, points=None, periods=1):
		if points == None:
			points = self.Ns
		points = float(points)
		times = n.arange(0., self.T * periods + self.T * periods / points, self.T * periods / points)
		return times
