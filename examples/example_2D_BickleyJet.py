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

sim = Simulation()  # Create a simulation object
sim.run_name = '2D_Bickley_Jet'

# Geometry and Model Equations
sim.geomx       = 'periodic'       # Geometry Types: 'periodic' or 'walls'
sim.geomy       = 'walls'
sim.stepper     = Step.AB3         # Time-stepping algorithm: Euler, AB2, RK4
sim.method      = 'Spectral'       # Numerical method: 'Spectral'
sim.dynamics    = 'Nonlinear'      # Dynamics: 'Nonlinear' or 'Linear'
sim.flux_method = Flux.spectral_sw # Flux method: spectral_sw is only option currently

# Specify paramters
sim.Lx  = 200e3           # Domain extent               (m)
sim.Ly  = 200e3           # Domain extent               (m)
sim.Nx  = 64             # Grid points in x
sim.Ny  = 64             # Grid points in y
sim.Nz  = 1               # Number of layers
sim.g   = 9.81            # Gravity                     (m/sec^2)
sim.f0  = 1.e-4           # Coriolis                    (1/sec)
sim.beta = 0e-10          # Coriolis beta               (1/m/sec)
sim.Hs  = [100.]          # Vector of mean layer depths (m)
sim.rho = [1025.]         # Vector of layer densities   (kg/m^3)
sim.end_time = 500*24.*hour   # End Time                    (sec)

# Parallel? Only applies to the FFTWs
sim.num_threads = 32

# Plotting parameters
sim.plott   = 12.*hour  # Period of plots
sim.animate = 'Save'      # 'Save' to create video frames,
                          # 'Anim' to animate,
                          # 'None' otherwise
sim.plot_vars = ['vort', 'v', 'u', 'h']
sim.clims = [ [-0.8, 0.8], [-0.5,0.5], [], []]                         

# Output parameters
sim.output = True        # True or False
sim.savet  = 10.*day      # Time between saves

# Diagnostics parameters
sim.diagt    = 2.*minute  # Time for output
sim.diagnose = False      # True or False

# Initialize the grid and zero solutions
sim.initialize()

for ii in range(sim.Nz):  # Set mean depths
    sim.soln.h[:,:,ii] = sim.Hs[ii]

# Bickley Jet initial conditions
# First we define the jet
Ljet = 10e3            # Jet width
amp  = 0.1             # Elevation of free-surface in basic state
sim.soln.h[:,:,0] += -amp*np.tanh(sim.Y/Ljet)
sim.soln.u[:,:,0]  =  sim.g*amp/(sim.f0*Ljet)/(np.cosh(sim.Y/Ljet)**2)

# Then we add on a random perturbation
np.random.seed(seed=100)
sim.soln.u[:,:,0] +=  1e-3*np.exp(-(sim.Y/Ljet)**2)*np.random.randn(sim.Nx,sim.Ny)

sim.run()                # Run the simulation


