'''
Created on Aug 5, 2013

@author: anne
'''
import numpy as n
class testlinbc:
	def leftBc(self,error):
		"""Velocity prescribed"""
		i=0
		rhoblk,Ublk,pblk=self.blocks(i)
		dUdt=1j*self.omega*self.U()[i]
		drhodt=1j*self.omega*self.rho()[i]
		dUdx=(self.U()[i+1]-self.U()[i])/(self.x[i+1]-self.x[i])
		error[rhoblk[0]:rhoblk[1]]=drhodt+1.2*dUdx
		error[rhoblk[0]]=self.rho()[i,0]-1.2
	#FIXTHIS
		error[Ublk[0]:Ublk[1]]=self.U()[i]
		error[Ublk[0]+1]=self.U()[i,1]-1. #Piston velocity
		error[Ublk[0]+2]=self.U()[i,2]-0.5 #Piston velocity
		#Implementation of equation of state
		error[pblk[0]:pblk[1]]=self.p()[i]-287*1.4*293.15*self.rho()[i]
		error[pblk[0]]=self.p()[i,0]-101325.
		return error
		return error

	def rightBc(self,error):
		#Velocity equals zero
		i=self.gridpoints-1
		rhoblk,Ublk,pblk=self.blocks(i)
		
		dUdt=1j*self.omega*self.U()[i]
		drhodt=1j*self.omega*self.rho()[i]
		dUdx=(self.U()[i]-self.U()[i-1])/(self.x[i]-self.x[i-1])
		error[rhoblk[0]:rhoblk[1]]=drhodt+1.2*dUdx
		error[rhoblk[0]]=self.rho()[i,0]-1.2
	#FIXTHIS
		error[Ublk[0]:Ublk[1]]=self.U()[i]
		#Implementation of equation of state
		error[pblk[0]:pblk[1]]=self.p()[i]-287*1.4*293.15*self.rho()[i]
		error[pblk[0]]=self.p()[i,0]-101325.
		return error