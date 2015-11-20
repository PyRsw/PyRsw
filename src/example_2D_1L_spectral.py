import numpy as np
import matplotlib.pyplot as plt

import Steppers as Step
import Fluxes as Flux
from PyRsw import Simulation
from constants import minute, hour, day
import sys

sim = Simulation()  # Create a simulation object

# Geometry and Model Equations
sim.geomx       = 'periodic'       # Geometry Types: 'periodic' or 'walls'
sim.geomy       = 'walls'
sim.stepper     = Step.AB2         # Time-stepping algorithm: Euler, AB2, RK4
sim.method      = 'Spectral'       # Numerical method: 'Spectral'
sim.dynamics    = 'Nonlinear'      # Dynamics: 'Nonlinear' or 'Linear'
sim.flux_method = Flux.spectral_sw # Flux method: spectral_sw is only option currently

# Specify paramters
sim.Lx  = 4000e3          # Domain extent               (m)
sim.Ly  = 4000e3          # Domain extent               (m)
sim.Nx  = 128             # Grid points in x
sim.Ny  = 128             # Grid points in y
sim.Nz  = 1               # Number of layers
sim.g   = 9.81            # Gravity                     (m/sec^2)
sim.f0  = 1.e-4           # Coriolis                    (1/sec)
sim.beta = 1e-11          # Coriolis beta               (1/m/sec)
sim.cfl = 0.01            # CFL coefficient             (m)
sim.Hs  = [100.]          # Vector of mean layer depths (m)
sim.rho = [1025.]         # Vector of layer densities   (kg/m^3)
sim.end_time = 4*24.*hour   # End Time                    (sec)

# Parallel? Only applies to the FFTWs
sim.num_threads = 4

# Plotting parameters
sim.plott   = 15.*minute  # Period of plots
sim.animate = 'Save'      # 'Save' to create video frames,
                          # 'Anim' to animate,
                          # 'None' otherwise
                         
# Output parameters
sim.output = False        # True or False
sim.savet  = 1.*hour      # Time between saves

# Diagnostics parameters
sim.diagt    = 2.*minute  # Time for output
sim.diagnose = False      # True or False

# Initialize the grid and zero solutions
sim.initialize()

for ii in range(sim.Nz):  # Set mean depths
    sim.soln.h[:,:,ii] = sim.Hs[ii]

# Gaussian initial conditions
x0 = 1.*sim.Lx/2.      # Centre
W  = 200.e3                # Width
amp = 0.1                  # Amplitude
sim.soln.h[:,:,0] += amp*np.exp(-(sim.x-x0)**2/(W**2)).reshape((sim.Nx,1))*np.exp(-(sim.y-x0)**2/(W**2)).reshape((1,sim.Ny))
#sim.soln.h[:,:,0] += amp*np.random.randn(sim.Nx,sim.Ny)

sim.run()                # Run the simulation


