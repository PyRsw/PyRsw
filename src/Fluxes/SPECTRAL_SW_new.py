import Differentiation as Diff
import numpy as np
import sys

def spectral_sw_flux(sim):

    # Difference in layer deformations gives layer depths
    hs = np.zeros((sim.Nx,sim.Ny,sim.Nz+1))
    hs[:,:,:sim.Nz] = sim.soln.h[:,:,:sim.Nz] - sim.soln.h[:,:,1:]

    # Loop through each layer and compute the flux 
    for ii in range(sim.Nz):

        h = hs[:,:,ii].reshape((sim.Nx,sim.Ny))
        u = sim.soln.u[:,:,ii].reshape((sim.Nx,sim.Ny))
        v = sim.soln.v[:,:,ii].reshape((sim.Nx,sim.Ny))

        # Coriolis terms
        sim.curr_flux.u[:,:,ii] = -sim.f0*v
        sim.curr_flux.v[:,:,ii] =  sim.f0*u

        # x fluxes
        if sim.Nx > 1:
            # Compute derivatives
            du = sim.ddx_u(u,sim.ik)
            dv = sim.ddx_v(v,sim.ik)
            dh = sim.ddx_h(h,sim.ik)

            # Intra-layer dyanics
            sim.curr_flux.u[:,:,ii] += -u*du - sim.gs[ii]*dh
            sim.curr_flux.v[:,:,ii] += -u*dv    
            sim.curr_flux.h[:,:,ii] += -h*du - u*dh

        # y fluxes
        if sim.Ny > 1:
            # Compute derivatives
            du = sim.ddy_u(u,sim.il)
            dv = sim.ddy_v(v,sim.il)
            dh = sim.ddy_h(h,sim.il)

            # Intra-layer dynamics
            sim.curr_flux.u[:,:,ii] += -v*du
            sim.curr_flux.v[:,:,ii] += -v*dv - sim.gs[ii]*dh
            sim.curr_flux.h[:,:,ii] += -h*dv - v*dh

    return

"""
def spectral_sw_lin_flux(sim):

    # Difference in layer deformations gives layer depths
    hs = np.zeros((sim.Nx,sim.Ny,sim.Nz+1))
    hs[:,:,:sim.Nz] = sim.soln.h[:,:,:sim.Nz] - sim.soln.h[:,:,1:]

    # Loop through each layer and compute the flux
    for ii in range(sim.Nz):

        h = hs[:,:,ii].reshape((sim.Nx,sim.Ny))
        u = sim.soln.u[:,:,ii].reshape((sim.Nx,sim.Ny))
        v = sim.soln.v[:,:,ii].reshape((sim.Nx,sim.Ny))

        # Coriolis terms
        sim.curr_flux.u[:,:,ii] = -sim.f0*v
        sim.curr_flux.v[:,:,ii] =  sim.f0*u

        # x fluxes
        if sim.Nx > 1:
            # Compute derivatives
            du = sim.ddx_u(u,sim.ik)
            dv = sim.ddx_v(v,sim.ik)
            dh = sim.ddx_h(h,sim.ik)

            # Intra-layer dyanics
            sim.curr_flux.u[:,:,ii] += -sim.gs[ii]*dh
            sim.curr_flux.h[:,:,ii] += -sim.Hs[ii]*du 

        # y fluxes
        if sim.Ny > 1:
            # Compute derivatives
            du = sim.ddy_u(u,sim.il)
            dv = sim.ddy_v(v,sim.il)
            dh = sim.ddy_h(h,sim.il)

            # Intra-layer dynamics
            sim.curr_flux.v[:,:,ii] += -sim.gs[ii]*dh
            sim.curr_flux.h[:,:,ii] += -sim.Hs[ii]*dv

    return
"""

def spectral_sw_source(sim):
    return 

def spectral_sw(sim):
    sim.x_derivs = Diff.SPECTRAL_x
    sim.y_derivs = Diff.SPECTRAL_y
    sim.flux_function = spectral_sw_flux
    sim.source_function = spectral_sw_source
