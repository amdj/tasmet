'''
Created on Aug 5, 2013

@author: anne
'''
import numpy as n
class testlinbc:
	def leftBc(self,error):
		"""Velocity prescribed"""
		i=0
		upiston=n.zeros((self.rho.shape[1]))
		upiston[1]=1.
		rhoblk,Ublk,pblk=self.blocks(i)
		V0=(self.x[i+1]-self.x[i])*self.Sf[0]
		dmdt=1j*self.omega*self.rho()[i]*V0+1.2*upiston*self.Sf[i]
		Udotnds=0.5*(self.Sf[i+1]*self.U()[i+1]+self.U()[i]*self.Sf[i])

		error[rhoblk[0]:rhoblk[1]]=dmdt+1.2*Udotnds
		error[rhoblk[0]]=self.rho()[i,0]-1.2
	#FIXTHIS
		dmfdt=1.2*1j*self.omega*self.U()[i]*V0
		pdotn=-self.Sf[i]*self.p()[i]+0.5*(self.Sf[i]*self.p()[i]+self.Sf[i+1]*self.p()[i+1])
		error[Ublk[0]:Ublk[1]]=dmfdt+pdotn
#		error[Ublk[0]+1]=dmfdt+pdotn #Piston velocity
#		error[Ublk[0]+2]=self.U()[i,2]-0.5 #Piston velocity
		#Implementation of equation of state
		c_sq=287*1.4*293.15
		error[pblk[0]:pblk[1]]=self.p()[i]-c_sq*self.rho()[i]
		error[pblk[0]]=self.p()[i,0]-101325.
		return error
	def rightBc2(self,error):
		#Velocity equals zero
		i=self.gridpoints-1
		rhoblk,Ublk,pblk=self.blocks(i)
		Vi=self.Sf[i]*(self.x[i]-self.x[i-1])
		Udotn=-0.5*(self.U()[i]+self.U()[i-1])
		drhodt=1j*self.omega*self.rho()[i]*Vi+1.2*Udotn

		dUdt=1j*self.omega*self.U()[i]
		error[rhoblk[0]:rhoblk[1]]=drhodt+1.2*Udotn
		error[rhoblk[0]]=self.rho()[i,0]-1.2

		pdotn=self.Sf[i]*self.p()[i]-0.5*(self.Sf[i]*self.p()[i]+self.Sf[i-1]*self.p()[i-1])
	#FIXTHIS
		error[Ublk[0]:Ublk[1]]=1.2*dUdt+pdotn
		#Implementation of equation of state
		c_sq=287*1.4*293.15
		error[pblk[0]:pblk[1]]=self.p()[i]-c_sq*self.rho()[i]
		error[pblk[0]]=self.p()[i,0]-101325.
		return error
	def rightBc(self,error):
		#Velocity equals zero
		i=self.gridpoints-1
		rhoblk,Ublk,pblk=self.blocks(i)


		drhodt=1j*self.omega*self.rho()[i]
		Udotn=0.5*(self.U()[i]+self.U()[i-1])
		error[rhoblk[0]:rhoblk[1]]=drhodt+1.2*Udotn
		error[rhoblk[0]]=self.rho()[i,0]-1.2
	#FIXTHIS
		pdotn=0.5*(self.p()[i]*self.Sf[i]+self.p()[i-1]*self.Sf[i-1])
		dUdt=1j*self.omega*self.U()[i]
		error[Ublk[0]:Ublk[1]]=1.2*dUdt+pdotn
		#Implementation of equation of state

		error[pblk[0]:pblk[1]]=self.p()[i]-287*1.4*293.15*self.rho()[i]
		error[pblk[0]]=self.p()[i,0]-101325.
		return error