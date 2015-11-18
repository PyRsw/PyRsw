import Differentiation as Diff
import numpy as np
import matplotlib.pyplot as plt
import sys

try:
    from pyfftw.interfaces.scipy_fftpack import fft, ifft, fftfreq, fftn, ifftn
except:
    from scipy.fftpack import fft, ifft, fftfreq, fftn, ifftn

def spectral_sw_flux(sim):

    # Loop through each layer and compute the flux
    for ii in range(sim.Nz):

        h = sim.soln.h[:,:,ii].reshape((sim.Nx,sim.Ny))
        u = sim.soln.u[:,:,ii].reshape((sim.Nx,sim.Ny))
        v = sim.soln.v[:,:,ii].reshape((sim.Nx,sim.Ny))

        # Compute x-derivatives
        du, dv, dh = sim.ddx_u(u,sim), sim.ddx_v(v,sim), sim.ddx_h(h,sim)

        # Coriolis and x-derivatives
        sim.curr_flux.u[:,:,ii] = - u*du + sim.f0*v  - sim.gs[ii]*dh
        sim.curr_flux.v[:,:,ii] = - u*dv - sim.f0*u     
        sim.curr_flux.h[:,:,ii] = - u*dh - h*du     

        # Compute y-derivatives
        du, dv, dh = sim.ddy_u(u,sim), sim.ddy_v(v,sim), sim.ddy_h(h,sim)
        
        # y-derivatives
        sim.curr_flux.u[:,:,ii] += - v*du
        sim.curr_flux.v[:,:,ii] += - v*dv - sim.gs[ii]*dh
        sim.curr_flux.h[:,:,ii] += - v*dh - h*dv 

        """
        print " u "
        print sim.curr_flux.u[:,:,ii] 
        print " v "
        print sim.curr_flux.v[:,:,ii] 
        print " h "
        print sim.curr_flux.h[:,:,ii] 
        sys.exit()
        """
        
    return

def spectral_sw_linear_flux(sim):

    # Loop through each layer and compute the flux
    for ii in range(sim.Nz):

        h = sim.soln.h[:,:,ii].reshape((sim.Nx,sim.Ny))
        u = sim.soln.u[:,:,ii].reshape((sim.Nx,sim.Ny))
        v = sim.soln.v[:,:,ii].reshape((sim.Nx,sim.Ny))

        # Compute x-derivatives
        du, dv, dh = sim.ddx_u(u,sim), sim.ddx_v(v,sim), sim.ddx_h(h,sim)

        # Coriolis and x-derivatives
        sim.curr_flux.u[:,:,ii] =   sim.f0*v - sim.gs[ii]*dh
        sim.curr_flux.v[:,:,ii] = - sim.f0*u     
        sim.curr_flux.h[:,:,ii] = - sim.Hs[ii]*du     

        # Compute y-derivatives
        du, dv, dh = sim.ddy_u(u,sim), sim.ddy_v(v,sim), sim.ddy_h(h,sim)
        
        # y-derivatives
        sim.curr_flux.v[:,:,ii] += - sim.gs[ii]*dh
        sim.curr_flux.h[:,:,ii] += - sim.Hs[ii]*dv

    return

def filter_general(sim):
    for ii in range(sim.Nz):
        ue = sim.soln.u[:,:,ii]
        ve = sim.soln.v[:,:,ii]
        he = sim.soln.h[:,:,ii]

        # Extend Grid if walls in x
        if sim.geomx=='walls':
            ue = np.concatenate([ue,-ue[::-1,:]],axis=0)
            ve = np.concatenate([ve, ve[::-1,:]],axis=0)
            he = np.concatenate([he, he[::-1,:]],axis=0)

        # Extend Grid if walls in y
        if sim.geomy=='walls':
            ue = np.concatenate([ue, ue[:,::-1]],axis=1)
            ve = np.concatenate([ve,-ve[:,::-1]],axis=1)
            he = np.concatenate([he, he[:,::-1]],axis=1)

        # Filter
        ue = ifftn(sim.sfilt*fftn(ue,axes=[0,1]),axes=[0,1]).real
        ve = ifftn(sim.sfilt*fftn(ve,axes=[0,1]),axes=[0,1]).real
        he = ifftn(sim.sfilt*fftn(he,axes=[0,1]),axes=[0,1]).real

        # Project on physical space
        sim.soln.u[:,:,ii] = ue[0:sim.Nx,0:sim.Ny]
        sim.soln.v[:,:,ii] = ve[0:sim.Nx,0:sim.Ny]
        sim.soln.h[:,:,ii] = he[0:sim.Nx,0:sim.Ny]
        
def spectral_sw(sim):

    if sim.Nx == 1:
        sim.Nkx = 1
    else:
        if sim.geomx == 'periodic':
            sim.Nkx = sim.Nx
        elif sim.geomx == 'walls':
            sim.Nkx = 2*sim.Nx
            
    if sim.Ny == 1:
        sim.Nky = 1
    else:
        if sim.geomy == 'periodic':
            sim.Nky = sim.Ny
        elif sim.geomy == 'walls':
            sim.Nky = 2*sim.Ny
        
    sim.x_derivs = Diff.SPECTRAL_x
    sim.y_derivs = Diff.SPECTRAL_y

    if sim.dynamics == 'Nonlinear':
        sim.flux_function = spectral_sw_flux
    elif sim.dynamics == 'Linear':
        sim.flux_function = spectral_sw_linear_flux
    else:
        print "dynamics must be from the list: Nonlinear, Linear"
        sys.exit()
            
    sim.apply_filter = filter_general

