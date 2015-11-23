import numpy as np
import matplotlib.pyplot as plt

import Steppers as Step
import Fluxes as Flux
from PyRsw import Simulation
from constants import minute, hour, day
import sys

sim = Simulation()  # Create a simulation object

# Geometry and Model Equations
sim.geomx       = 'walls'       # Geometry Types: 'periodic' or 'walls'
sim.geomy       = 'periodic'
sim.stepper     = Step.AB3         # Time-stepping algorithm: Euler, AB2, RK4
sim.method      = 'Spectral'       # Numerical method: 'Spectral'
sim.dynamics    = 'Nonlinear'      # Dynamics: 'Nonlinear' or 'Linear'
sim.flux_method = Flux.spectral_sw # Flux method: spectral_sw is only option currently

# Specify paramters
sim.Lx  = 4000e3          # Domain extent               (m)
sim.Ly  = 4000e3          # Domain extent               (m)
sim.Nx  = 1                 # Grid points in x
sim.Ny  = 128             # Grid points in y
sim.Nz  = 1               # Number of layers
sim.g   = 9.81            # Gravity                     (m/sec^2)
sim.f0  = 1.e-4           # Coriolis                    (1/sec)
sim.beta = 0e-10          # Coriolis beta parameter     (1/m/sec)
sim.cfl = 0.2             # CFL coefficient             (m)
sim.Hs  = [100.]          # Vector of mean layer depths (m)
sim.rho = [1025.]         # Vector of layer densities   (kg/m^3)
sim.end_time = 14*24.*hour   # End Time                    (sec)

# Plotting parameters
sim.plott   = 20.*minute  # Period of plots
sim.animate = 'Anim'      # 'Save' to create video frames,
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
W  = 200.e3                # Width
amp = 1.                  # Amplitude
if sim.Ny==1:
    sim.soln.h[:,:,0] += amp*np.exp(-(sim.X)**2/(W**2))
elif sim.Nx==1:
    sim.soln.h[:,:,0] += amp*np.exp(-(sim.Y)**2/(W**2))

sim.run()                # Run the simulation


# Hovmuller plot
plt.figure()
t = np.arange(0,sim.end_time+sim.plott,sim.plott)/86400.

if sim.Ny==1:
    x = sim.x/1e3
elif sim.Nx == 1:
    x = sim.y/1e3

for L in range(sim.Nz):
    field = sim.hov_h[:,0,:].T - np.sum(sim.Hs[L:])
    cv = np.max(np.abs(field.ravel()))
    plt.subplot(sim.Nz,1,L+1)
    plt.pcolormesh(x,t, field,
        cmap=sim.cmap, vmin = -cv, vmax = cv)
    plt.axis('tight')
    plt.title(r"$\mathrm{Hovm{\"o}ller} \; \mathrm{Plot} \; \mathrm{of} \; \eta$", fontsize = 16)
    if sim.Nx > 1:
        plt.xlabel(r"$\mathrm{x} \; \mathrm{(km)}$", fontsize=14)
    else:
        plt.xlabel(r"$\mathrm{y} \; \mathrm{(km)}$", fontsize=14)
    plt.ylabel(r"$\mathrm{Time} \; \mathrm{(days)}$", fontsize=14)
    plt.colorbar()

plt.show()
