import Differentiation as Diff
import numpy as np
import sys

def spectral_sw_xy_flux(sim):

    # Difference in layer deformations gives layer depths
    hs = np.zeros((sim.Nx,sim.Ny,sim.Nz+1))
    hs[:,:,:sim.Nz] = sim.soln.h[:,:,:sim.Nz] - sim.soln.h[:,:,1:]

    # Loop through each layer and compute the flux
    for ii in range(sim.Nz):

        h = hs[:,:,ii].reshape((sim.Nx,sim.Ny))
        u = sim.soln.u[:,:,ii].reshape((sim.Nx,sim.Ny))
        v = sim.soln.v[:,:,ii].reshape((sim.Nx,sim.Ny))

        # Compute x-derivatives
        du, dv, dh = sim.ddx_u(u,sim.ik), sim.ddx_v(v,sim.ik), sim.ddx_h(h,sim.ik)

        # Coriolis and x-derivatives
        sim.curr_flux.u[:,:,ii] = -sim.f0*v - u*du - sim.gs[ii]*dh
        sim.curr_flux.v[:,:,ii] =  sim.f0*u -u*dv    
        sim.curr_flux.h[:,:,ii] = -h*du - u*dh

        # Compute y-derivatives
        du, dv, dh = sim.ddy_u(u,sim.il), sim.ddy_v(v,sim.il), sim.ddy_h(h,sim.il)
        
        # y-derivatives
        sim.curr_flux.u[:,:,ii] += -v*du
        sim.curr_flux.v[:,:,ii] += -v*dv - sim.gs[ii]*dh
        sim.curr_flux.h[:,:,ii] += -h*dv - v*dh

    return

def spectral_sw_lin_xy_flux(sim):

    # Difference in layer deformations gives layer depths
    hs = np.zeros((sim.Nx,sim.Ny,sim.Nz+1))
    hs[:,:,:sim.Nz] = sim.soln.h[:,:,:sim.Nz] - sim.soln.h[:,:,1:]

    # Loop through each layer and compute the flux
    for ii in range(sim.Nz):

        h = hs[:,:,ii].reshape((sim.Nx,sim.Ny))
        u = sim.soln.u[:,:,ii].reshape((sim.Nx,sim.Ny))
        v = sim.soln.v[:,:,ii].reshape((sim.Nx,sim.Ny))

        # Compute x-derivatives
        du, dv, dh = sim.ddx_u(u,sim.ik), sim.ddx_v(v,sim.ik), sim.ddx_h(h,sim.ik)

        # Coriolis and x-derivatives
        sim.curr_flux.u[:,:,ii] = -sim.f0*v - sim.gs[ii]*dh
        sim.curr_flux.v[:,:,ii] =  sim.f0*u 
        sim.curr_flux.h[:,:,ii] = -sim.Hs[ii]*du 

        # Compute y-derivatives
        du, dv, dh = sim.ddy_u(u,sim.il), sim.ddy_v(v,sim.il), sim.ddy_h(h,sim.il)
        
        # y-derivatives
        sim.curr_flux.v[:,:,ii] += -sim.gs[ii]*dh
        sim.curr_flux.h[:,:,ii] += -sim.Hs[ii]*dv

    return

def spectral_sw_x_flux(sim):

    # Difference in layer deformations gives layer depths
    hs = np.zeros((sim.Nx,sim.Ny,sim.Nz+1))
    hs[:,:,:sim.Nz] = sim.soln.h[:,:,:sim.Nz] - sim.soln.h[:,:,1:]

    # Loop through each layer and compute the flux
    for ii in range(sim.Nz):

        h = hs[:,:,ii].reshape((sim.Nx,sim.Ny))
        u = sim.soln.u[:,:,ii].reshape((sim.Nx,sim.Ny))
        v = sim.soln.v[:,:,ii].reshape((sim.Nx,sim.Ny))

        # Compute x-derivatives
        du, dv, dh = sim.ddx_u(u,sim.ik), sim.ddx_v(v,sim.ik), sim.ddx_h(h,sim.ik)

        # Coriolis and x-derivatives
        sim.curr_flux.u[:,:,ii] = -sim.f0*v - u*du - sim.gs[ii]*dh
        sim.curr_flux.v[:,:,ii] =  sim.f0*u -u*dv    
        sim.curr_flux.h[:,:,ii] = -h*du - u*dh

    return

def spectral_sw_lin_x_flux(sim):

    # Difference in layer deformations gives layer depths
    hs = np.zeros((sim.Nx,sim.Ny,sim.Nz+1))
    hs[:,:,:sim.Nz] = sim.soln.h[:,:,:sim.Nz] - sim.soln.h[:,:,1:]

    # Loop through each layer and compute the flux
    for ii in range(sim.Nz):

        h = hs[:,:,ii].reshape((sim.Nx,sim.Ny))
        u = sim.soln.u[:,:,ii].reshape((sim.Nx,sim.Ny))
        v = sim.soln.v[:,:,ii].reshape((sim.Nx,sim.Ny))

        # Compute x-derivatives
        du, dv, dh = sim.ddx_u(u,sim.ik), sim.ddx_v(v,sim.ik), sim.ddx_h(h,sim.ik)

        # Coriolis and x-derivatives
        sim.curr_flux.u[:,:,ii] = -sim.f0*v - sim.gs[ii]*dh
        sim.curr_flux.v[:,:,ii] =  sim.f0*u     
        sim.curr_flux.h[:,:,ii] = -sim.Hs[ii]*du 

    return


def spectral_sw_y_flux(sim):

    # Difference in layer deformations gives layer depths
    hs = np.zeros((sim.Nx,sim.Ny,sim.Nz+1))
    hs[:,:,:sim.Nz] = sim.soln.h[:,:,:sim.Nz] - sim.soln.h[:,:,1:]

    # Loop through each layer and compute the flux
    for ii in range(sim.Nz):

        h = hs[:,:,ii].reshape((sim.Nx,sim.Ny))
        u = sim.soln.u[:,:,ii].reshape((sim.Nx,sim.Ny))
        v = sim.soln.v[:,:,ii].reshape((sim.Nx,sim.Ny))

        # Compute y-derivatives
        du, dv, dh = sim.ddy_u(u,sim.il), sim.ddy_v(v,sim.il), sim.ddy_h(h,sim.il)
        
        # y-derivatives
        sim.curr_flux.u[:,:,ii] = - sim.f0*v - v*du
        sim.curr_flux.v[:,:,ii] =   sim.f0*u - v*dv - sim.gs[ii]*dh
        sim.curr_flux.h[:,:,ii] = - h*dv - v*dh

    return

def spectral_sw_lin_y_flux(sim):

    # Difference in layer deformations gives layer depths
    hs = np.zeros((sim.Nx,sim.Ny,sim.Nz+1))
    hs[:,:,:sim.Nz] = sim.soln.h[:,:,:sim.Nz] - sim.soln.h[:,:,1:]

    # Loop through each layer and compute the flux
    for ii in range(sim.Nz):

        h = hs[:,:,ii].reshape((sim.Nx,sim.Ny))
        u = sim.soln.u[:,:,ii].reshape((sim.Nx,sim.Ny))
        v = sim.soln.v[:,:,ii].reshape((sim.Nx,sim.Ny))

        # Compute y-derivatives
        du, dv, dh = sim.ddy_u(u,sim.il), sim.ddy_v(v,sim.il), sim.ddy_h(h,sim.il)
        
        # y-derivatives
        sim.curr_flux.u[:,:,ii] = - sim.f0*v 
        sim.curr_flux.v[:,:,ii] =   sim.f0*u - sim.gs[ii]*dh
        sim.curr_flux.h[:,:,ii] = - sim.Hs[ii]*dv

    return

def spectral_sw(sim):
    sim.x_derivs = Diff.SPECTRAL_x
    sim.y_derivs = Diff.SPECTRAL_y

    if sim.geom == 'x':
        if sim.dynamics == 'Nonlinear':
            sim.flux_function = spectral_sw_x_flux
        if sim.dynamics == 'Linear':
            sim.flux_function = spectral_sw_lin_x_flux
            
    if sim.geom == 'y':
        if sim.dynamics == 'Nonlinear':
            sim.flux_function = spectral_sw_y_flux
        if sim.dynamics == 'Linear':
            sim.flux_function = spectral_sw_lin_y_flux
            
    if sim.geom == 'xy':
        if sim.dynamics == 'Nonlinear':
            sim.flux_function = spectral_sw_xy_flux
        if sim.dynamics == 'Linear':
            sim.flux_function = spectral_sw_lin_xy_flux

