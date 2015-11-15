import Steppers as Step
import Fluxes as Flux
import numpy as np
from PyRsw import Simulation
from constants import minute, hour, day

# Create a simulation object
sim = Simulation()

# Specify goemetry: 'x', 'y', 'xy'
sim.geom  = 'x' 

# Specify geometry conditions: 'periodic' or 'wall'
sim.geomx = 'periodic'

# Specify time-stepping algorithm: Euler, AB2, RK4
sim.stepper = Step.AB2

# Specify numerical method: 'Spectral' or 'FV'
sim.method = 'Spectral'

#FJP: add this somewhere
# Nonlinear or Linear
sim.dynamics = 'Linear'

#FJP: if want linear option better to have e instead of h
# Specify the flux method: spectral_sw is only option currently
sim.flux_method = Flux.spectral_sw

# Specify paramters
sim.Lx  = 2000e3          # Domain extent (m)
#sim.Ly  = 2000e3          # Domain extent (m)
sim.Nx  = 128*2             # Grid points in x
sim.Ny  = 1               # Grid points in y
sim.Nz  = 1               # Number of layers
sim.f0  = 4.e-4           # Coriolis
sim.cfl = 0.1             # CFL coefficient
sim.Hs  = [100.]          # Vector of mean layer depths
sim.rho = [1025.]         # Vector of layer densities
sim.end_time = 36.*hour   # End Time

#FJP: Anim + Hov
#FJP: NL vs Linear: see flux fcn
#FJP: change readthedocs to pyrsw not pysw

# Plotting parameters
sim.plott = 20.*minute   # Period of plots
sim.animate = 'Anim'     # 'Save' to create video frames,
                         # 'Anim' to animate,
                         # 'None' otherwise
#FJP: push this later
sim.ylims[2] = [-1,1]    # Manual ylimits on plots: ylims[0] -> u, ylims[1] -> v, ylim[2] -> h
                         #    Specifies the ylimits for line plots and the climits for 2d pcolor

# FJP: ???
# Output parameters
sim.output = False      # True or False
sim.savet  = 1.*hour    # Time between saves

# FJP: ???
# Diagnostics parameters
sim.diagt = 2.*minute
sim.diagnose = False

# Initialize the grid and zero solutions
sim.initialize()

# Set mean depths
for ii in range(sim.Nz):
    sim.soln.h[:,:,ii] = np.sum(sim.Hs[ii:])

#print sim.Nx, sim.Ny
#print len(sim.x), len(sim.y)

# Specify initial conditions to be Gaussian
x0 = 1.*sim.Lx/2.      # Centre
W  = 50.e3             # Width
amp = 1.               # Amplitue
sim.soln.h[:,:,0] += amp*np.exp(-(sim.x-x0)**2/(W**2)).reshape((sim.Nx,1))
#sim.soln.h[:,:,0] += amp*np.exp(-(sim.y-x0)**2/(W**2)).reshape((1,sim.Ny))

# Run the simulation
sim.run()
