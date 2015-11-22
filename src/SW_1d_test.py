#!/usr/bin/env python
#  SW_1d.m
#
# Solve the 1-Layer and 1D Rotating Shallow Water (SW) Model
#
# Fields: 
#   u : zonal velocity
#   v : meridional velocity
#   h : fluid depth
#
# Evolution Eqns:
#   B = g*h + 0.5*(u**2 + v**2)     Bernoulli function
#   Z = v_x + f                     Total Vorticity
#   q = Z/h                         Potential Vorticity
#   [U,V] = h[u,v]                  Transport velocities
#
#	u_t =  (q*V^x) + d_x h
#	v_t = -(q*U^y)^x 
#	h_t = - d_x[U]
#
# Geometry: periodic in x 
#           Non-staggered grids
#           All fields are at half grid points
#
# Numerical Method:
# 1) Fourier Spectral methods for derivatives
# 2) Exponential Filter in space
# 3) Adams-Bashforth for time stepping
#

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys

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

def plot_uvh(parms,uvh,t):

    index = parms.fplot
    
    plt.plot(parms.x/1e3,uvh[index,:])
    #plt.xlim([0,parms.Lx/1e3])
    #plt.ylim([0,hmax])
    name = parms.pname[index] + " at t = %5.2f days" % (t/3600./24.)
    plt.title(name)

#######################################################
#         Flux for SW                                 #
#######################################################

def ddx_period(f,ik): 

    df = np.real(ifft(ik*fft(f)))

    return df

def ddx_even(f,ik):

    N = len(f)
    fe = np.hstack((f,np.flipud(f)))
    df = np.real(ifft(ik*fft(fe)))
    df = df[0:N]
    
    return df

def ddx_odd(f,ik): 

    N = len(f)
    fe = np.hstack((f,-np.flipud(f)))
    df = np.real(ifft(ik*fft(fe)))
    df = df[0:N]
    
    return df

def flux_sw_periodic(uvh, parms):

    # Define parameters
    dx = parms.dx
    gp = parms.gp
    f0 = parms.f0 
    H0 = parms.H0
    nx = parms.nx
    ik = parms.ik
    
    # Compute Fields
    u = uvh[0,:]
    v = uvh[1,:]
    h = uvh[2,:] + H0
    B = gp*h + 0.5*(u**2 + v**2)

    #Derivative of Bernoulli function
    B_x = ddx_period(B,ik)

    # Potential Vorticity
    q = (ddx_period(v,ik) + f0)/h

    # Derivative of mass flux
    U = h*u
    U_x = ddx_period(U,ik)
    
    # Compute fluxes
    flux = np.vstack([ q*h*v - B_x, -q*U, -U_x])

    #compute energy and enstrophy
    energy = 0.5*np.mean( gp*h**2 + h*(u**2 + v**2) )
    enstrophy = 0.5*np.mean(q**2*h)

    return flux, energy, enstrophy

def flux_sw_walls(uvh, parms):

    # Define parameters
    dx = parms.dx
    gp = parms.gp
    f0 = parms.f0 
    H0 = parms.H0
    nx = parms.nx
    ik = parms.ik
    r  = 0*parms.r
    
    # Compute Fields
    u = uvh[0,:]
    v = uvh[1,:]
    h = uvh[2,:] + H0

    #Compute x-derivatives
    du = ddx_odd( u,ik)
    dv = ddx_even(v,ik)
    dh = ddx_even(h,ik)
    
    # Compute fluxes
    flux = np.vstack([ -u*du + f0*v - gp*dh,
                       -u*dv - f0*u,
                       -u*dh - h*du])

    #compute energy and enstrophy
    energy = 0 #0.5*np.mean( gp*h**2 + h*(u**2 + v**2) )
    enstrophy = 0 #0.5*np.mean(q**2*h)

    return flux, energy, enstrophy

def flux_sw_qB_walls(uvh, parms):

    # Define parameters
    dx = parms.dx
    gp = parms.gp
    f0 = parms.f0 
    H0 = parms.H0
    nx = parms.nx
    ik = parms.ik
    r  = 0*parms.r
    
    # Compute Fields
    u = uvh[0,:]
    v = uvh[1,:]
    h = uvh[2,:] + H0
    B1 = 0.5*u**2 + gp*h
    B2 = 0.5*v**2

    #Derivative of Bernoulli function
    B_x = ddx_even(B1,ik) + ddx_odd(B2,ik)

    # Potential Vorticity
    q = (ddx_even(v,ik) + f0)/h

    # Derivative of mass flux
    U = h*u
    U_x = ddx_odd(U,ik)
    
    # Compute fluxes
    flux = np.vstack([ q*h*v - B_x -r*u, -q*U, -U_x - r*uvh[2,:]])

    #compute energy and enstrophy
    energy = 0.5*np.mean( gp*h**2 + h*(u**2 + v**2) )
    enstrophy = 0.5*np.mean(q**2*h)

    return flux, energy, enstrophy

#######################################################
#        Parameters Class                             #
#######################################################

class Parms(object):
    """A class to solve the one-layer Sw model."""
    
    def __init__(
        self,
        # grid size parameters
        nx=128,                     # grid resolution
        Lx=4000e3,                  # zonal domain size 
        #geom='periodic',            # type of geometry: periodic, walls
        geom='walls',               # type of geometry: periodic, walls
                 
        # physical parameters
        gp  = 9.81,                 # (reduced) gravity
        H0  = 1000,                 # mean depth
        f0  = 0e-4,                 # Coriolis parameter
        beta= 0e-11,                # gradient of coriolis parameter
        r0  = 0.0,                  # Coefficient for sponge layer

        # timestepping parameters
        dt=60.,                     # numerical timstep
        tplot=600.,                 # interval for plots (in timesteps)
        tmax=10*86400.,             # total time of integration

        # what field to plot
        fplot=2,
        ):
        
        # put all the parameters into the object grid
        self.nx = nx
        self.Lx = Lx
        self.geom = geom

        if geom =='periodic':
            method = flux_sw_periodic

            # Compute wavenumbers
            kx = 2*np.pi/Lx*np.hstack([range(0,int(nx/2)), range(-int(nx/2),0)])
        elif geom == 'walls':
            method = flux_sw_walls
            
            # Compute wavenumbers
            kx = np.pi/Lx*np.hstack([range(0,int(nx)), range(-int(nx),0)])
        else:
            print "Geometry must be from the list: periodic, walls"
            sys.exit()
        self.method = method        
        self.ik = 1j*kx

        # physical parameters
        self.gp   = gp
        self.H0   = H0
        self.f0   = f0
        self.beta = beta
        self.r0   = r0

        # timestepping
        self.dt = dt
        self.tplot = tplot
        self.npt = int(tplot/dt)
        self.tmax = tmax
        self.nt = int(tmax/dt)
        
        # Define Grid (staggered C-grid)
        dx = Lx/nx
        self.dx = dx
        self.x = np.linspace(0+dx/2,Lx-dx/2,nx)

        self.pname = ['u','v','eta']
        self.fplot = fplot

        # Define sponge
        r = 0*self.x
        spI = nx/8
        delta = self.x[-1] - self.x[-spI]
        r[-spI:] = ((self.x[-spI:] - self.x[-spI])/delta)**2
        r[0:(spI-1)] = ((self.x[0:(spI-1)] - self.x[spI-1])/delta)**2
        self.r = np.sqrt(gp*H0)/dx*r
        
        # Filter Parameters
        #kmax = max(kx);
        #ks = 0.4*kmax;  
        #km = 0.5*kmax;
        #alpha = 0.69*ks**(-1.88/np.log(km/ks));
        #beta  = 1.88/np.log(km/ks);
        #self.sfilt = np.exp(-alpha*(kx**2)**(beta/2.0));
        
        fcut, ford, fstr = 0.6, 2.0, 20.0
        k = kx/max(kx.ravel())
        self.sfilt = np.exp(-fstr*((np.abs(k)-fcut)/(1-fcut))**ford)*(np.abs(k)>fcut) + (np.abs(k)<fcut)
        #self.sfilt = filtx.reshape((len(k),1))
        

#######################################################
#        Solve SW model                               #
#######################################################

def solve_sw(parms, uvh0):
    
    # initialize fields
    uvh       = np.zeros(parms.nx)
    energy    = np.zeros(parms.nt)
    enstrophy = np.zeros(parms.nt)    
    
    # Euler Step
    t,ii = 0., 0
    NLm, energy[ii], enstrophy[ii] = parms.method(uvh0, parms)
    uvh = uvh0 + parms.dt*NLm;

    for ii in range(2,parms.nt+1):

        # AB2 step
        t = (ii-1)*parms.dt
        NL, energy[ii], enstrophy[ii] = parms.method(uvh, parms)
        uvh = uvh + 0.5*parms.dt*(3*NL - NLm)

        # Extend Grid if walls in x
        if (parms.geom=='walls'):
            ue = np.concatenate([uvh[0,:],-uvh[0,::-1]],axis=0)
            ve = np.concatenate([uvh[1,:], uvh[1,::-1]],axis=0)
            he = np.concatenate([uvh[2,:], uvh[2,::-1]],axis=0)

        # Filter
        ue = ifft(parms.sfilt*fft(ue)).real
        ve = ifft(parms.sfilt*fft(ve)).real
        he = ifft(parms.sfilt*fft(he)).real
        
        # Project on physical space
        uvh[0,:] = ue[0:parms.nx]
        uvh[1,:] = ve[0:parms.nx]
        uvh[2,:] = he[0:parms.nx]

        # Reset fluxes
        NLm = NL

        if (ii-0)%parms.npt==0:

            print "t  = %f hours norm is %f" % (t/3600./24., np.amax(uvh[2,:]))

            plt.clf()
            plot_uvh(parms,uvh,t)
            #plt.ylim([-0.2,1.2*hmax])
            #plt.show(block=False)
            plt.pause(0.01)
            plt.draw()
            #plt.show()
            #sys.exit()

    return uvh, energy, enstrophy

def solve_sw_AB3(parms, uvh0):
    
    # initialize fields
    uvh       = np.zeros(parms.nx)
    energy    = np.zeros(parms.nt)
    enstrophy = np.zeros(parms.nt)    
    
    # Euler Step
    t,ii = 0., 0
    NLnm, energy[ii], enstrophy[ii] = parms.method(uvh0, parms)
    uvh = uvh0 + parms.dt*NLnm;

    # AB2 step
    t,ii = parms.dt, 1
    NLn, energy[ii], enstrophy[ii] = parms.method(uvh, parms)
    uvh = uvh + 0.5*parms.dt*(3*NLn - NLnm)

    for ii in range(3,parms.nt+1):

        #FJP: change filter???
        # AB3 step
        t = (ii-1)*parms.dt
        NL, energy[ii-1], enstrophy[ii-1] = parms.method(uvh, parms);
        uvh  = uvh + parms.dt/12*(23*NL - 16*NLn + 5*NLnm)

        # Extend Grid if walls in x
        if (parms.geom=='walls'):
            ue = np.concatenate([uvh[0,:],-uvh[0,::-1]],axis=0)
            ve = np.concatenate([uvh[1,:], uvh[1,::-1]],axis=0)
            he = np.concatenate([uvh[2,:], uvh[2,::-1]],axis=0)

        # Filter
        ue = ifft(parms.sfilt*fft(ue)).real
        ve = ifft(parms.sfilt*fft(ve)).real
        he = ifft(parms.sfilt*fft(he)).real
        
        # Project on physical space
        uvh[0,:] = ue[0:parms.nx]
        uvh[1,:] = ve[0:parms.nx]
        uvh[2,:] = he[0:parms.nx]

        # Reset fluxes
        NLnm = NLn
        NLn  = NL

        if (ii-0)%parms.npt==0:

            print "t  = %f hours norm is %f" % (t/3600./24., np.amax(uvh[2,:]))

            plt.clf()
            plot_uvh(parms,uvh,t)
            #plt.ylim([-0.2,1.2*hmax])
            #plt.show(block=False)
            plt.pause(0.01)
            plt.draw()
            #plt.show()
            #sys.exit()

    return uvh, energy, enstrophy


#######################################################
#         Main Program                                #
#######################################################

#FJP: AB3 can set dt = 10/sc for numerical stability
#FJP: AB2 can set dt =  5/sc for numerical stability 
# Numerical parameters
sc = 1
nx = 128*sc
Lx = 200e3
dt = 2./sc
tmax = 1.*3600.*24.
geom = 'walls'
fplot = 2

# Physical parameters
f0, gp, H0  = 1e-4, 9.81, 100

# Set parameters
parms = Parms(f0=f0, gp=gp, H0=H0, Lx=Lx, nx=nx, dt=dt, tmax=tmax, geom=geom, fplot=fplot)
nx = parms.nx

# Jet initial conditions
Lj = 10.e3                 # Width
amp = 0.1                  # Amplitude
uvh0 = np.zeros((3,nx))
uvh0[1,:]  =  parms.gp*amp/(parms.f0*Lj)/(np.cosh((parms.x-parms.Lx/2)/Lj))**2
uvh0[2,:]  =  amp*np.tanh((parms.x-parms.Lx/2)/Lj)
uvh0[2,:] += 1e-5*np.random.randn(parms.nx)

# Plot Initial Conditions
plt.ion()
plot_uvh(parms,uvh0,0)
#plt.xlim([0,parms.Lx])
#plt.ylim([-0.2,1.2*hmax])
#plt.show(block=False)
plt.pause(0.01)
plt.draw()
#plt.show()

q,energy,enstrophy = solve_sw(parms, uvh0)

plt.ioff()
plt.show()

sys.exit()

print "Error in energy is ", np.amax(energy-energy[0])
print "Error in enstrophy is ", np.amax(enstrophy-enstrophy[0])
    
fig, axarr = plt.subplots(2, sharex=True)
ax1 = plt.subplot(2,1,1)
ax2 = plt.subplot(2,1,2)
ax1.plot((energy-energy[0]),'-ob',linewidth=2, label='Energy')
ax1.set_title('Energy')
ax2.plot((enstrophy-enstrophy[0]),'-or', linewidth=2, label='Enstrophy')
ax2.set_title('Enstrophy')
plt.show()    
            
sys.exit()

