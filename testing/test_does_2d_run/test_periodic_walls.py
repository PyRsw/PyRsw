import numpy as np
import matplotlib.pyplot as plt
import sys

# Add the PyRsw tools to the path
# At the moment it is given explicitely.
# In the future, it could also be added to the
# pythonpath environment variable
sys.path.append('../src')

import Steppers as Step
import Fluxes as Flux
from PyRsw import Simulation
from constants import minute, hour, day

def test():

    sim = Simulation()  # Create a simulation object

    # Geometry and Model Equations
    sim.geomx       = 'periodic'
    sim.geomy       = 'walls'
    sim.stepper     = Step.AB3      
    sim.method      = 'Spectral'    
    sim.dynamics    = 'Nonlinear'  
    sim.flux_method = Flux.spectral_sw 

    # Specify paramters
    sim.Lx  = 4000e3  
    sim.Ly  = 4000e3   
    sim.Nx  = 128       
    sim.Ny  = 128      
    sim.Nz  = 1         
    sim.g   = 9.81       
    sim.f0  = 1.e-4       
    sim.beta = 0e-11       
    sim.cfl = 0.1           
    sim.Hs  = [100.]         
    sim.rho = [1025.]         
    sim.end_time = 5.*minute 

    sim.animate = 'None'    
    sim.output = False       
    sim.diagnose = False

    # Initialize the grid and zero solutions
    sim.initialize()

    for ii in range(sim.Nz):  # Set mean depths
        sim.soln.h[:,:,ii] = sim.Hs[ii]

    # Gaussian initial conditions
    W  = 200.e3                # Width
    amp = 1.                  # Amplitude
    sim.soln.h[:,:,0] += amp*np.exp(-(sim.X/W)**2 - (sim.Y/W)**2)

    # Run the simulation
    sim.run() 

