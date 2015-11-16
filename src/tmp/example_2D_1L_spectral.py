import Steppers as Step
import Fluxes as Flux
import numpy as np
from PySW import Simulation
from constants import minute, hour, day

# Create a simulation object
sim = Simulation()

# Specify numerical method: 'Spectral' or 'FV'
sim.method = 'Spectral'

# Specify geometry conditions: 'periodic' or 'wall'
sim.geomx = 'periodic'
sim.geomy = 'periodic'

# Specify time-stepping algorithm: Euler, AB2, RK4
sim.stepper = Step.AB2

# Specify the flux method: spectral_sw is only option currently 
sim.flux_method = Flux.spectral_sw

# Specify paramters
sim.Lx  = 1000e3     # Domain extent (m)
sim.Ly  = 600e3     # Domain extent (m)
sim.Nx  = 128*4     # Grid points
sim.Ny  = 128*4
sim.Nz  = 1         # Number of layers
sim.f0  = 1.e-4     # Coriolis
sim.cfl = 0.2      # CFL coefficient
sim.Hs  = [50.]     # Vector of mean layer depths
sim.rho = [1025.]   # Vector of layer densities

# Plotting parameters
sim.plott = 8.*minute   # Time between plots
sim.animate = 'Save'    # 'Save' to create video frames, 'Anim' to animate, 'None' otherwise
sim.ylims[2] = [48,52]  # Manual ylimits on plots: ylims[0] -> u, ylims[1] -> v, ylim[2] -> h

# Output parameters
sim.output = False      # True or False
sim.savet  = 1.*hour    # Time between saves

# Diagnostics parameters
sim.diagt = 2.*minute
sim.diagnose = False

# Initialize the grid and zero solutions
sim.initialize()

# Specify initial conditions
for ii in range(sim.Nz):
    sim.soln.h[:,:,ii] = np.sum(sim.Hs[ii:])

x0 = 1.*sim.Lx/2.
y0 = 1.*sim.Ly/2.
W  = 50.e3 
amp = 1. 
X,Y = np.meshgrid(sim.x,sim.y,indexing='ij')
sim.soln.h[:,:,0] += amp*np.exp(-((X-x0)/W)**2 - ((Y-y0)/W)**2).reshape((sim.Nx,sim.Ny))

# Run the simulation
sim.end_time = 32.*hour
sim.run()
