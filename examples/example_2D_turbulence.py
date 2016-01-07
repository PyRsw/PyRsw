import numpy as np
import matplotlib.pyplot as plt
import sys

from scipy.fftpack import fft, ifft, fftfreq, fftn, ifftn, fftshift

# Add the PyRsw tools to the path
# At the moment it is given explicitely.
# In the future, it could also be added to the
# pythonpath environment variable
sys.path.append('../src')

import Steppers as Step
#import Fluxes as Flux
from PyRsw import Simulation
from constants import minute, hour, day

sim = Simulation()  # Create a simulation object
sim.run_name = '2D_turbulence_N64'

# Geometry and Model Equations
sim.geomx       = 'periodic'       # Geometry Types: 'periodic' or 'walls'
sim.geomy       = 'periodic'
sim.stepper     = Step.AB3         # Time-stepping algorithm: Euler, AB2, RK4
sim.method      = 'Sadourny'       # Numerical method: 'Spectral'
sim.dynamics    = 'Nonlinear'      # Dynamics: 'Nonlinear' or 'Linear'
#sim.flux_method = Flux.SADOURNY_SW

# Specify paramters
sim.Lx  = 400e3           # Domain extent               (m)
sim.Ly  = 400e3           # Domain extent               (m)
sim.Nx  = 64             # Grid points in x
sim.Ny  = 64             # Grid points in y
sim.Nz  = 1               # Number of layers
sim.g   = 5e-2            # Gravity                     (m/sec^2)
sim.f0  = 1.e-4           # Coriolis                    (1/sec)
sim.beta = 0e-10          # Coriolis beta               (1/m/sec)
sim.Hs  = [4000.]          # Vector of mean layer depths (m)
sim.rho = [1025.]         # Vector of layer densities   (kg/m^3)
sim.end_time = 2000*24.*hour   # End Time                    (sec)

# Parallel? Only applies to the FFTWs
sim.num_threads = 32

# Plotting parameters
sim.plott   = 5*24*hour  # Period of plots
sim.animate = 'Save'      # 'Save' to create video frames,
                          # 'Anim' to animate,
                          # 'None' otherwise
sim.plot_vars = ['vort']
sim.clims = [ [-.1, .1] ]                         
#sim.plot_vars = ['vort', 'v', 'u', 'h']
#sim.clims = [ [-0.8, 0.8], [-0.5,0.5], [], []]                         

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

# Random initial conditions
#
# Spectral Grid
kx = 2*np.pi/sim.Lx*np.hstack([range(0,int(sim.Nx/2)), range(-int(sim.Nx/2),0)])
ky = 2*np.pi/sim.Ly*np.hstack([range(0,int(sim.Ny/2)), range(-int(sim.Ny/2),0)])
KX, KY = np.meshgrid(kx, ky, indexing='ij')
K = np.sqrt(KX**2 + KY**2)

K0 = kx[10]
W = 2*(kx[1] - kx[0])
amp  = 200.0             # Elevation of free-surface in basic state

theta = 2.*np.pi*np.random.randn(sim.Nx,sim.Ny)
etah = 1/K*np.exp(-((K-K0)/W)**2)*np.exp(1j*theta)
etah[0,0] = 0.0
eta = ifftn(etah).real
eta = eta/np.max(eta)*amp
sim.soln.h[:-1,:-1,0] += eta
sim.soln.h[-1,:,0] = sim.soln.h[0,:,0]
sim.soln.h[:,-1,0] = sim.soln.h[:,0,0]
#sim.soln.h[:,:,0] += eta

#plt.clf()
#plt.contourf(fftshift(KX),fftshift(KY), fftshift(np.abs(etah)))
#plt.pcolormesh(sim.grid_x.h, sim.grid_y.h, sim.soln.h[:,:,0])
#plt.colorbar()
#plt.show()
#sys.exit()

sim.run()                # Run the simulation


