#!/usr/bin/python
#
# Import the TASMET models in a 'ta' namespace
import TASMET as ta

# Create some global configuration 
Nf=6 # Number of harmonics to solve for
freq=100 # Fundamental oscillating frequency [Hz]
gas='air' # The working gas
gc=ta.Globalconf(Nf,freq,gas)

# Create segments, such as tubes, segments, etc
seg1=...
seg2=...

# Create connectors and boundary conditions to put
# all stuff together
con1=...
con2=...

# Create a TaSystem object
sys=ta.TaSystem(gc)

# Next step: add the segments and connectors to the system
sys+=seg1
sys+=seg2
sys+=...
sys+=con1
sys+=con2

# Next step: solve the system
sol=ta.Solver()
# Sometimes, the solver requires some extra parameters, such
# as convergence criteria, and the maximum number of iterations
sc=ta.SolverConfiguration(param1=val1,param2=val2)

# Set the solverconfiguration to the solver:
sol.setSc(sc)

############### The big step: solve the system! ############
sol.Solve(sys)

############### If the solver is done, do postprocessing.

# Postprocessing code here...
