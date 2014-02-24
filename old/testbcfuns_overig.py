'''
Created on Aug 5, 2013

@author: anne
'''
	def leftBc1(self,error):
		i=0
		rhoblk,Ublk,pblk,Tblk=self.blocks(i)
		#Vi=self.Sf[i]*2*self.x[i]
		#dmdt=1j*self.omega*self.rho()[i]*Vi
		error[Ublk[0]:Ublk[1],0]=self.U()[i]-self.testU
		error[pblk[0]:pblk[1],0]=self.p()[i]-self.testp
		error[Tblk[0]:Tblk[1],0]=self.T()[i]-self.testT
		error[rhoblk[0]:rhoblk[1],0]=self.rho()[i]-self.dft(self.m.rho(self.p.getTdata()[i],self.T.getTdata()[i]))
		return error	
	def leftBc2(self,error):
		i=0
		rhoblk,Ublk,pblk,Tblk=self.blocks(i)
		Vi=self.Sf[i]*2*self.x[i]

		Uleft_td=self.idft1long(self.testU)
		uleft_td=Uleft_td/self.SfT[i]
		minleft_td=Uleft_td*self.rho.getTdata()[i]
		minleft=self.dft(minleft_td)
		rhoudotn=0.5*(self.mf()[i]+self.mf()[i+1])-minleft
		dmdt=1j*self.omega*self.rho()[i]*Vi
		#Density equation
		error[rhoblk[0]:rhoblk[1],0]=dmdt+rhoudotn
		error[rhoblk[0],0]=self.rho.getAdata()[0,0]-1.2

		#Momentum balance
		dmfdt=1j*self.omega*self.mf()[i]*Vi
		mominleft_td=uleft_td*minleft_td
		mominleft=self.dft(mominleft_td)
		#print mominleft
		rhouudotn=0.5*(self.momf()[i]+self.momf()[i+1])-mominleft
		pdotn=-self.Sf[i]*self.p()[i]+0.5*(self.Sf[i]*self.p()[i]+self.Sf[i+1]*self.p()[i+1])
		error[Ublk[0]:Ublk[1],0]=dmfdt+pdotn+rhouudotn

		dEtotdt=1j*self.omega*self.Etot()[i]*Vi
		Hleft_td=minleft_td*(self.m.h_T(self.T.getTdata()[i])+0.5*uleft_td*uleft_td)
		Hleft=self.dft(Hleft_td)

		error[Tblk[0]:Tblk[1],0]=dEtotdt+0.5*(self.Hf()[i+1]-self.Hf()[i])-Hleft
		error[Tblk[0],0]=self.T()[i,0]-293.15
		error[pblk[0]:pblk[1],0]=self.p()[i]-self.dft(self.m.p_rho_T(self.rho.getTdata()[i], self.T.getTdata()[i]))

		return error