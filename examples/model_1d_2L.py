#!/usr/bin/env python
#  pysw_1d_2L.py
#
# Solve the 1D, 2-LayerRotating Shallow Water (SW) Model
# FJP: make multi-layer
# FJP: make 2D
#
# Fields: 
#   u : zonal velocity
#   v : meridional velocity
#   h : depth of each layer
#   e : free-surface deformations
#
# Evolution Eqns:
#	u_t = - v*u_y + f*v
#	v_t = - v*v_y - f*u - h_y 
#	h_t = - (h*v)_y
#
# Geometry: y: periodic or slip walls 
#           Non-staggered grids
#           All fields are at half grid points
#
# Numerical Method:
# 1) Fourier Spectral methods for derivatives
# 2) Exponential Filter in space
# 3) Adams-Bashforth for time stepping
# 4) Uses a "pretty good sponge"
# 5) Forces a surface gravity wave

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys

# FJP: import our own FFTW wrappers with f2py?
# https://sysbio.ioc.ee/projects/f2py2e/usersguide.html#tth_sEcB.1
# FJP: or use Cython like in pyqg

# Put sponge back in?
# Linear vs Nonlinear option?
# Look at AMATH900 for details on power spectrum
# Make Hovemoeller plots

try:
    import pyfftw
    from numpy import zeros as nzeros

    # Keep fft objects in cache for efficiency
    nthreads = 1
    pyfftw.interfaces.cache.enable()
    pyfftw.interfaces.cache.set_keepalive_time(1e8)
    def empty(N, dtype="float", bytes=16):
        return pyfftw.n_byte_align_empty(N, bytes, dtype=dtype)

    def zeros(N, dtype="float", bytes=16):
        return pyfftw.n_byte_align(nzeros(N, dtype=dtype), bytes)
    
    # Monkey patches for fft
    ifft = pyfftw.interfaces.numpy_fft.ifft
    fft = pyfftw.interfaces.numpy_fft.fft

except:    
    from scipy.fftpack import fft, ifft
    print Warning("Install pyfftw, it is much faster than numpy fft")

# Allows us to store parameters in one variable
class parameters:
    pass

#######################################################
#        Plotting Functions                           #
#######################################################

# FJP: needs some polishing
def plot_uvh(parms,uvh,t):

    index = parms.fplot #2 - corrosponds to eta
    
    plt.plot(parms.y/1e3,uvh[index,:])
    plt.xlim([0,parms.Ly/1e3])
    #plt.ylim([0,hmax])
    plt.ylim([-1.5e-2, 1.5e-2])
    name = parms.pname[index] + " at t = %5.2f days" % (t/3600./24.)
    plt.title(name)
    

#######################################################
#         Flux for SW                                 #
#######################################################

# FJP: have simple switch between methods
# spectra, sadourney, weno5?
# periodic or Neumann bcs?
# diffusion or not

def ddy_period(f,ik):

    df = np.real(ifft(ik*fft(f)))

    return df

def ddy_even(f,ik):

    N = len(f)
    fe = np.hstack((f,np.flipud(f)))
    df = np.real(ifft(ik*fft(fe)))
    df = df[0:N]
    
    return df

def ddy_odd(f,ik):

    N = len(f)
    fe = np.hstack((f,-np.flipud(f)))
    df = np.real(ifft(ik*fft(fe)))
    df = df[0:N]
    
    return df

def flux_sw_periodic(uvh, parms, t):

    # Define parameters
    dy = parms.dy
    gp = parms.gp
    f0 = parms.f0 
    H0 = parms.H0
    ny = parms.ny
    f  = parms.f
    ik = parms.ik
    
    # Compute Fields
    u = uvh[0,:]
    v = uvh[1,:]
    h = uvh[2,:] + H0

    # Compute derivatives
    u_y = ddy_period(u,ik)
    v_y = ddy_period(v,ik)
    h_y = ddy_period(h,ik)

    # Compute fluxes
    flux = np.vstack([ -v*u_y + f*v,
                       -v*v_y - f*u - gp*h_y - sp*v,
                       -h*v_y - v*h_y - sp*uvh[2,:]])

    #compute energy and enstrophy
    # FJP: correct
    energy = 0.5*np.mean( gp*h**2 + h*(u**2 + v**2) )
    enstrophy = 0.5*np.mean( (-u_y + f)**2/h)

    return flux, energy, enstrophy

# FJP: combine this with above
# FJP: remove tidal forcing
def flux_sw_walls(uvh, parms, t):

    # Define parameters
    dy = parms.dy
    gp = parms.gp
    f0 = parms.f0 
    H0 = parms.H0
    ny = parms.ny
    f  = parms.f
    ik = parms.ik
        
    # Compute Fields
    u = uvh[0,:]
    v = uvh[1,:]
    h = uvh[2,:] + H0
    
    # Compute derivatives
    u_y = ddy_even(u,ik)
    v_y = ddy_odd(v,ik)
    h_y = ddy_even(h,ik)

    # Compute fluxes
    flux = np.vstack([ -v*u_y + f*v, 
                       -v*v_y - f*u - h_y,
                       -h*v_y - v*h_y])
    
    #compute energy and enstrophy
    energy = 0.5*np.mean( gp*h**2 + h*(u**2 + v**2) )
    enstrophy = 0.5*np.mean( (-u_y + f)**2/h)

    return flux, energy, enstrophy

#######################################################
#        Parameters Class                             #
#######################################################

class Parms(object):
    """A class to solve the one-layer Sw model."""

    # Default parameters
    def __init__(
        self,
        # grid size parameters
        ny=128,                     # grid resolution
        Ly=4000e3,                  # zonal domain size 
        geom='walls',               # type of geometry: periodic, walls
                 
        # physical parameters
        gp  = 9.81,                 # (reduced) gravity
        H0  = 1000,                 # mean depth
        f0  = 0e-4,                 # Coriolis parameter
        beta= 0e-11,                # gradient of coriolis parameter

        # timestepping parameters
        dt=60.,                     # numerical timstep
        tplot=600.,              # interval for plots (in timesteps)
        tmax=10*86400.,             # total time of integration
    ):
        
        # put all the parameters into the object grid
        self.ny = ny
        self.Ly = Ly
        self.geom = geom

        if geom =='periodic':
            method = flux_sw_periodic

            # Compute wavenumbers
            ky = 2*np.pi/Ly*np.hstack([range(0,int(ny/2)), range(-int(ny/2),0)])
        elif geom == 'walls':
            method = flux_sw_walls
            
            # Compute wavenumbers
            ky = np.pi/Ly*np.hstack([range(0,int(ny)), range(-int(ny),0)])
        else:
            print "Geometry must be from the list: periodic, walls"
            sys.exit()
        self.method = method        
        self.ik = 1j*ky

        # physical parameters
        self.gp   = gp
        self.H0   = H0
        self.f0   = f0
        self.beta = beta

        # timestepping
        self.dt = dt
        self.tplot = tplot
        self.npt = int(tplot/dt)
        self.tmax = tmax
        self.nt = int(tmax/dt)
        
        # Define Grid (staggered C-grid)
        dy = Ly/ny
        self.dy = dy
        self.y = np.linspace(0+dy/2,Ly-dy/2,ny)

        # Coriolis parameter
        self.f = f0 + beta*self.y.copy()
        
        # Plotting parameters
        self.pname = ['u','v','eta']
        self.fplot = 2

        # Filter Parameters
        kmax = max(ky);
        ks = 0.4*kmax;
        km = 0.5*kmax;
        alpha = 0.69*ks**(-1.88/np.log(km/ks));
        beta  = 1.88/np.log(km/ks);
        self.sfilt = np.exp(-alpha*(ky**2)**(beta/2.0));


#######################################################
#        Solve SW model                               #
#######################################################

def solve_sw(parms, uvh0):
    
    # initialize fields
    uvh       = np.zeros(parms.ny)
    energy    = np.zeros(parms.nt)
    enstrophy = np.zeros(parms.nt)    

    # Build matrix to save
    hsave = np.zeros((parms.tmax/parms.tplot+1,parms.ny))
    cnt = 0
    hsave[cnt,:] = uvh0[2,:]
     
    # Euler Step
    t,ii = 0., 0
    NLnm, energy[ii], enstrophy[ii] = parms.method(uvh0, parms, t)
    uvh = uvh0 + parms.dt*NLnm;

    # AB2 step
    t,ii = parms.dt, 1
    NLn, energy[ii], enstrophy[ii] = parms.method(uvh, parms, t)
    uvh = uvh + 0.5*parms.dt*(3*NLn - NLnm)

    print t, np.linalg.norm(uvh)
    
    for ii in range(3,parms.nt+1):

        # AB3 step
        t = (ii-1)*parms.dt
        NL, energy[ii-1], enstrophy[ii-1] = parms.method(uvh, parms, t);
        uvh  = uvh + parms.dt/12*(23*NL - 16*NLn + 5*NLnm)

        #print t, np.linalg.norm(uvh)
        
        # Reset fluxes
        NLnm = NLn
        NLn  = NL
        
        if (ii-0)%parms.npt==0:

            print "t  = %f hours" % (t/3600./24.)

            hsave[cnt,:] = uvh[2,:]
            cnt += 1
            
            plt.clf()
            plot_uvh(parms,uvh,t)
            plt.show() #(block=False)
            plt.draw()
            #plt.show()
            #sys.exit()

    return hsave, energy, enstrophy

#######################################################
#         Main Program                                #
#######################################################

# Numerical parameters
sc = 1
ny = 128*sc
dt = 60./sc
tmax = 2*3600.*24.
geom = 'walls'

# Physical parameters
f0, gp, H0  = 1e-4, 9.81, 1000

# Set parameters
parms = Parms(f0=f0, gp=gp, H0=H0, ny=ny, dt=dt, tmax=tmax, geom=geom)
ny = parms.ny

print "ny = ", ny

# Specify initial conditions
t = 0
uvh0 = np.zeros((3,ny))
hmax = 1e-2
uvh0[2,:] = hmax*np.exp(-((parms.y-parms.Ly/2)**2)/(parms.Ly/20)**2)

# Plot Initial Conditions
plt.ion()
plot_uvh(parms,uvh0,0)
#plt.xlim([0,parms.Ly])
##plt.ylim([0,hmax])
#plt.show() #(block=False)
#plt.draw()
plt.show()
#sys.exit()

q,energy,enstrophy = solve_sw(parms, uvh0)

plt.ioff()
plt.show()

plt.clf()
plt.pcolormesh(tpt/(3600),parms.y/1e3,q.T)
plt.ylabel('y (km)')
plt.xlabel('Time (hours)')
plt.title('Forcing: 1D 1-layer SW')
plt.xlim([24, 48])
plt.colorbar()
plt.clim([-0.004, 0.004])
plt.show()

#sys.exit()

#print "Error in energy is ", np.amax(energy-energy[0])
#print "Error in enstrophy is ", np.amax(enstrophy-enstrophy[0])
    
#fig, axarr = plt.subplots(2, sharex=True)
#ax1 = plt.subplot(2,1,1)
#ax2 = plt.subplot(2,1,2)
#ax1.plot((energy-energy[0]),'-ob',linewidth=2, label='Energy')
#ax1.set_title('Energy')
#ax2.plot((enstrophy-enstrophy[0]),'-or', linewidth=2, label='Enstrophy')
#ax2.set_title('Enstrophy')
#plt.show()    
            

