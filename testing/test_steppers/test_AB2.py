import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('../src')

import Steppers as Step
import Fluxes as Flux
from PyRsw import Simulation
from constants import minute, hour, day


def test():

    sim = Simulation() 

    # Geometry and Model Equations
    sim.geomy       = 'periodic'  
    sim.stepper     = Step.AB2     
    sim.method      = 'Spectral'     
    sim.dynamics    = 'Nonlinear'     
    sim.flux_method = Flux.spectral_sw 

    # Specify paramters
    sim.Ly  = 4000e3   
    sim.Ny  = 128       
    sim.cfl = 0.5        
    sim.Hs  = [100.]      
    sim.rho = [1025.]      
    sim.end_time = 5.*minute

    # Plotting parameters
    sim.animate = 'None'  
    sim.output = False   
    sim.diagnose = False 

    # Initialize the grid and zero solutions
    sim.initialize()

    for ii in range(sim.Nz): 
        sim.soln.h[:,:,ii] = sim.Hs[ii]

    # Gaussian initial conditions
    x0 = 1.*sim.Lx/2.     
    W  = 200.e3          
    amp = 1.            
    sim.soln.h[:,:,0] += amp*np.exp(-(sim.Y)**2/(W**2))

    sim.run()      
   
