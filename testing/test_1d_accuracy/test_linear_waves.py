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
    sim.stepper     = Step.AB3      
    sim.method      = 'Spectral'     
    sim.dynamics    = 'Linear'     
    sim.flux_method = Flux.spectral_sw 

    # Specify paramters
    sim.Ly  = 4000e3   
    sim.Ny  = 256       
    sim.f0  = 0.
    sim.Hs  = [100.]      
    sim.rho = [1025.]      
    sim.end_time = sim.Ly/(np.sqrt(sim.Hs[0]*sim.g))

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
    IC = sim.soln.h[:,:,0].copy()

    sim.run()       

    # Compare final state to initial conditions
    error_h = np.linalg.norm(IC - sim.soln.h[:,:,0])
    error_v = np.linalg.norm(sim.soln.v[:,:,0])
    assert (error_h < 2e-5) and (error_v < 1e-7)

test()
