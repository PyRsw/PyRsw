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
sim.Lx  = 200e3           # Domain extent               (m)
sim.Ly  = 200e3           # Domain extent               (m)
sim.Nx  = 256             # Grid points in x
sim.Ny  = 256             # Grid points in y
sim.Nz  = 1               # Number of layers
sim.g   = 9.81            # Gravity                     (m/sec^2)
sim.f0  = 1.e-4           # Coriolis                    (1/sec)
sim.beta = 0e-10          # Coriolis beta               (1/m/sec)
sim.cfl = 0.2             # CFL coefficient             (m)
sim.Hs  = [100.]          # Vector of mean layer depths (m)
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
Ljet = 10e3            # Jet width
amp  = 0.1             # Elevation of free-surface in basic state
sim.soln.h[:,:,0] += -amp*np.tanh(sim.Y/Ljet)
sim.soln.u[:,:,0]  =  sim.g*amp/(sim.f0*Ljet)/(np.cosh(sim.Y/Ljet)**2)
sim.soln.u[:,:,0] +=  1e-3*np.exp(-(sim.Y/Ljet)**2)*np.random.randn(sim.Nx,sim.Ny)




# Define cheb array
def cheb(N):
    if N == 0:
        D = 0
        x = 1
    else:
        x = np.cos(np.pi*np.array(range(0,N+1))/N).reshape([N+1,1])
        c = np.ravel(np.vstack([2, np.ones([N-1,1]), 2])) \
            *(-1)**np.ravel(np.array(range(0,N+1)))
        c = c.reshape(c.shape[0],1)
        X = np.tile(x,(1,N+1))
        dX = X-(X.conj().transpose())
        D  = (c*(1/c).conj().transpose())/(dX+(np.eye(N+1)))   # off-diagonal entries
        D  = D - np.diag(np.sum(D,1))   # diagonal entries
    return D,x

## Specify wavenumber
Nk  = 200
kks = np.linspace(0,0.2*np.pi/sim.Ly,Nk+1)[1:-1]
k = kks[50]

# Define Differentiation Matrix and grid
Dy,y = cheb(sim.Ny)
y   = (y[:,0]+1.)*sim.Ly/2
Dy = Dy*(2/sim.Ly)
 
# Define Basic State
UB  = sim.g*amp/(sim.f0*Ljet)/(np.cosh(y/Ljet)**2)
dUB = np.dot(Dy,UB) 
HB  = sim.Hs[0] - amp*np.tanh(y/Ljet)

# Build matrices
Z1 = np.zeros((sim.Ny+1,sim.Ny+1))
I1 = np.eye(sim.Ny+1)

A = np.vstack(( np.hstack((          np.diag(UB,0), (np.diag(dUB - sim.f0))[:,1:-1],  sim.g*I1)),
                np.hstack(( -sim.f0/k**2*I1[1:-1,:],        np.diag(UB,0)[1:-1,1:-1], -sim.g/k**2*Dy[1:-1,:])),
                np.hstack((             np.diag(HB),                 (Dy*HB)[:,1:-1],  np.diag(UB,0)))))

# Compute linear stability
eigVals, eigVecs = spalg.eig(A)

# Sort eigenvalues and eigenvectors
ind = (-np.imag(eigVals)).argsort()
eigVecs = eigVecs[:,ind]
eigVals = k*eigVals[ind]

print np.real(eigVals[0:10])#*k*3600.*24.
print np.imag(eigVals[0:10])#*k*3600.*24.
#self.eigVecs[:,0:Ne,cnt] = eigVecs[:,0:Ne]
#self.eigVals[0:Ne,cnt]   = eigVals[0:Ne]

sys.exit()

#sim.run()                # Run the simulation


