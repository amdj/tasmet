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
    def __init__(self,n_particles,sys,steps=50,nperiods=50):
        self.particles=[]
        self.L=sys.getx()[-1]
        self.freq=sys.getFreq()
        for i in range(n_particles):
            self.particles.append(Particle(i/(n_particles-1)*self.L,self.freq,steps,nperiods))

    def getpos(self):
        pos=np.zeros(len(self.particles))
        for i,particle in enumerate(self.particles):
            pos[i]=particle.pos[0]
        return pos
            
    def u(t,x):
        return sys.quantityAtTimeAndPlaceInterp(t,x,"velo")
