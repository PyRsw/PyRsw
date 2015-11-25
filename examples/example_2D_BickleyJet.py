import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as spalg
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

# Geometry and Model Equations
sim.geomx       = 'periodic'       # Geometry Types: 'periodic' or 'walls'
sim.geomy       = 'walls'
sim.stepper     = Step.AB3         # Time-stepping algorithm: Euler, AB2, RK4
sim.method      = 'Spectral'       # Numerical method: 'Spectral'
sim.dynamics    = 'Nonlinear'      # Dynamics: 'Nonlinear' or 'Linear'
sim.flux_method = Flux.spectral_sw # Flux method: spectral_sw is only option currently

# Specify paramters
sim.Lx  = 350e3           # Domain extent               (m)
sim.Ly  = 350e3           # Domain extent               (m)
sim.Nx  = 256             # Grid points in x
sim.Ny  = 256             # Grid points in y
sim.Nz  = 1               # Number of layers
sim.g   = 9.81            # Gravity                     (m/sec^2)
sim.f0  = 1.e-4           # Coriolis                    (1/sec)
sim.beta = 0e-10          # Coriolis beta               (1/m/sec)
sim.cfl = 0.2             # CFL coefficient             (m)
sim.Hs  = [500.]          # Vector of mean layer depths (m)
sim.rho = [1025.]         # Vector of layer densities   (kg/m^3)
sim.end_time = 14*24.*hour   # End Time                    (sec)

# Parallel? Only applies to the FFTWs
sim.num_threads = 4

# Plotting parameters
sim.plott   = 30.*minute  # Period of plots
sim.animate = 'Save'      # 'Save' to create video frames,
                          # 'Anim' to animate,
                          # 'None' otherwise
sim.plot_vars = ['vort']
                         
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
Ljet = 20e3            # Jet width
Umax = 0.5
amp  = Umax*sim.f0*Ljet/sim.g    # Elevation of free-surface in basic state
sim.soln.h[:,:,0] += -amp*np.tanh(sim.Y/Ljet)
sim.soln.u[:,:,0]  =  sim.g*amp/(sim.f0*Ljet)/(np.cosh(sim.Y/Ljet)**2)
sim.soln.u[:,:,0] +=  1e-3*np.exp(-(sim.Y/Ljet)**2)*np.random.randn(sim.Nx,sim.Ny)

## Specify wavenumber
dkk = 2e-2
kk = np.arange(dkk,2+dkk,dkk)/Ljet
Nk = len(kk)

## Storage Arrays
# FJP: store eigenvectors too
ne = 2                           # number of eigenvalues to store
growsw = np.zeros((ne,Nk))
freqsw = np.zeros((ne,Nk))

# Define Differentiation Matrix and grid
from Stability_tools import cheb
Dy,y = cheb(sim.Ny)
y   = (y[:,0]+1)*sim.Ly/2
Dy = Dy*(2/sim.Ly)

# Define Basic State
sim.UB  =  sim.g*amp/(sim.f0*Ljet)/(np.cosh((y-sim.Ly/2.)/Ljet)**2)
sim.HB  = sim.Hs[0] - amp*np.tanh((y-sim.Ly/2.)/Ljet)
sim.dUB = np.dot(Dy,sim.UB) 

# Test Geostrophy
if np.amax(sim.f0*sim.UB+sim.g*np.dot(Dy,sim.HB)) > 1e-8:
    print "Geostrophic balance not satisfied.  Exit"
    sys.exit()    

# Compute spectrum
from Stability_tools import stability_sw
growsw, freqsw = stability_sw(sim, Dy, y, kk, ne)

print growsw*3600.*24

plt.clf
plt.plot(kk, 3600*24.*growsw[0,:],'ob')
plt.plot(kk, 3600*24.*growsw[1,:],'or')
plt.show()
sys.exit()

sim.run()                # Run the simulation


