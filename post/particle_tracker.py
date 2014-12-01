#!/usr/bin/python
import numpy as np
from scipy.integrate import ode

class Particle(object):
    def __init__(self,initpos,freq,steps=50,nperiods=5):
        self.initpos=initpos
        self.nperiods=nperiods
        self.T=1/freq
        self.dt=self.T/(steps)
        self.pos=np.zeros(steps*nperiods)
        self.t=np.zeros(steps*nperiods)
        self.pos[:]=initpos
    def integrate(self,ufun):
        r=ode(ufun).set_integrator('dopri5')
        r.set_initial_value(self.pos[0],0)
        i=1
        while r.successful() and i<self.t.size:
            r.integrate(r.t+self.dt)
            self.pos[i]=r.y
            self.t[i]=r.t
            i+=1
            
        
class ParticleTracker(object):
    def __init__(self,n_particles,sys,steps=50,nperiods=1):
        self.particles=[]
        self.sys=sys
        self.L=sys.getx()[-1]
        self.freq=sys.getFreq()
        dx=self.L/(n_particles)
        x=dx/2.
        for i in range(n_particles):
            self.particles.append(Particle(x,self.freq,steps,nperiods))
            print("Integrating particle %g..." %i)
            self.particles[i].integrate(self.u)
            x+=dx
        self.t=self.particles[0].t
    def getpos(self,ti=0):
        pos=np.zeros(len(self.particles))
        for i,particle in enumerate(self.particles):
            pos[i]=particle.pos[ti]
        return pos
            
    def u(self,t,x):
        return self.sys.quantityAtTimeAndPlaceInterp(t,x,"velo")
