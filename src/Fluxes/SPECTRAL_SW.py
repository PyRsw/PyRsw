import Differentiation as Diff
import numpy as np
import sys

try:
    import pyfftw
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
        sim.curr_flux.u[:,:,ii] = - u*du + sim.F*v  - sim.gs[ii]*dh
        sim.curr_flux.v[:,:,ii] = - u*dv - sim.F*u     
        sim.curr_flux.h[:,:,ii] = - u*dh - h*du     

        # Compute y-derivatives
        du, dv, dh = sim.ddy_u(u,sim), sim.ddy_v(v,sim), sim.ddy_h(h,sim)
        
        # y-derivatives
        sim.curr_flux.u[:,:,ii] += - v*du
        sim.curr_flux.v[:,:,ii] += - v*dv - sim.gs[ii]*dh
        sim.curr_flux.h[:,:,ii] += - v*dh - h*dv 

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
        sim.curr_flux.u[:,:,ii] =   sim.F*v - sim.gs[ii]*dh
        sim.curr_flux.v[:,:,ii] = - sim.F*u     
        sim.curr_flux.h[:,:,ii] = - sim.Hs[ii]*du     

        # Compute y-derivatives
        du, dv, dh = sim.ddy_u(u,sim), sim.ddy_v(v,sim), sim.ddy_h(h,sim)
        
        # y-derivatives
        sim.curr_flux.v[:,:,ii] += - sim.gs[ii]*dh
        sim.curr_flux.h[:,:,ii] += - sim.Hs[ii]*dv

    return

def filter_general(sim):
    for ii in range(sim.Nz):
        ue = sim.soln.u[:,:,ii].copy()
        ve = sim.soln.v[:,:,ii].copy()
        he = sim.soln.h[:,:,ii].copy()

        # Extend Grid if walls in x
        if (sim.geomx=='walls') and (sim.Nx > 1):
            ue = np.concatenate([ue,-ue[::-1,:]],axis=0)
            ve = np.concatenate([ve, ve[::-1,:]],axis=0)
            he = np.concatenate([he, he[::-1,:]],axis=0)

        # Extend Grid if walls in y
        if (sim.geomy=='walls') and (sim.Ny > 1):
            ue = np.concatenate([ue, ue[:,::-1]],axis=1)
            ve = np.concatenate([ve,-ve[:,::-1]],axis=1)
            he = np.concatenate([he, he[:,::-1]],axis=1)

        # Filter
        ue = sim.ifftn_u(sim.sfilt*sim.fftn_u(ue)).real
        ve = sim.ifftn_v(sim.sfilt*sim.fftn_v(ve)).real
        he = sim.ifftn_h(sim.sfilt*sim.fftn_h(he)).real

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

    # If possible, use pyfftw to preserve wisdom.
    try:
        tmp_in_u  = pyfftw.n_byte_align_empty((sim.Nkx,sim.Nky),16,dtype='complex128')
        tmp_out_u = pyfftw.n_byte_align_empty((sim.Nkx,sim.Nky),16,dtype='complex128')
        sim.fftn_u  = pyfftw.FFTW(tmp_in_u, tmp_out_u, axes=[0,1], direction='FFTW_FORWARD', threads = sim.num_threads)
        sim.ifftn_u = pyfftw.FFTW(tmp_out_u, tmp_in_u, axes=[0,1], direction='FFTW_BACKWARD', threads = sim.num_threads)

        tmp_in_v  = pyfftw.n_byte_align_empty((sim.Nkx,sim.Nky),16,dtype='complex128')
        tmp_out_v = pyfftw.n_byte_align_empty((sim.Nkx,sim.Nky),16,dtype='complex128')
        sim.fftn_v  = pyfftw.FFTW(tmp_in_v, tmp_out_v, axes=[0,1], direction='FFTW_FORWARD', threads = sim.num_threads)
        sim.ifftn_v = pyfftw.FFTW(tmp_out_v, tmp_in_v, axes=[0,1], direction='FFTW_BACKWARD', threads = sim.num_threads)

        tmp_in_h  = pyfftw.n_byte_align_empty((sim.Nkx,sim.Nky),16,dtype='complex128')
        tmp_out_h = pyfftw.n_byte_align_empty((sim.Nkx,sim.Nky),16,dtype='complex128')
        sim.fftn_h  = pyfftw.FFTW(tmp_in_h, tmp_out_h, axes=[0,1], direction='FFTW_FORWARD', threads = sim.num_threads)
        sim.ifftn_h = pyfftw.FFTW(tmp_out_h, tmp_in_h, axes=[0,1], direction='FFTW_BACKWARD', threads = sim.num_threads)

    except:
        def ifftx(ar):
            return ifftn(ar,axes=[0,1])
        def fftx(ar):
            return fftn(ar,axes=[0,1])

        sim.fftn_u  = fftn
        sim.ifftn_u = ifftn
        sim.fftn_v  = fftn
        sim.ifftn_v = ifftn
        sim.fftn_h  = fftn
        sim.ifftn_h = ifftn
        print('Failed to import pyfftw for filter transforms.')

        
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

