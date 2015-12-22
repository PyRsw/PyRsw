import Differentiation as Diff
import numpy as np
import sys

# Geometry: periodic in x and y
#           Arakawa C-grid
#
#      |           |          |         |
#      h --  u --  h  -- u -- h -- u -- h --
#      |           |          |         |
#      |           |          |         |
#      v     q     v     q    v    q    v
#      |           |          |         |
#      |           |          |         |
#      h --  u --  h  -- u -- h -- u -- h --
#      |           |          |         |
#      |           |          |         |
#      v     q     v     q    v    q    v
#      |           |          |         |
#      |           |          |         |
#      h --  u --  h  -- u -- h -- u -- h --
#      |           |          |         | 
#      |           |          |         |
#      v     q     v     q    v    q    v
#      |           |          |         |
#      |           |          |         |
#      h --  u --  h  -- u -- h -- u -- h --
#

# If Periodic X Periodic: u,v,h = Nx by Ny
# If Periodic X Walls:  u,h: Nx by (Ny+1) and
#                         v: Nx by Ny
# If Walls X Periodic:  v,h: (Nx+1) by Ny and
#                         u: Nx by Ny      
#
# N,S rows:
# -> must advance u,h
# -> Must extend v to compute V_y
# -> If v = 0 then maybe (q*V^x) = 0 too?

# W,E columns:
# -> must advance v,h
# -> Must extend u to compute U_x
# -> If u = 0 then (q*U^y) = 0 too?

# ghost cells:
# u-eqn: need q*V^x
# v-eqn: need q*U_y
# h-eqn: need U left and V down
def sadourny_sw_flux(sim):

    Nx, Ny, Nz = sim.Nx, sim.Ny, sim.Nz
    dx, dy     = sim.dx[0], sim.dx[1]
    
    # Loop through each layer and compute the flux
    for ii in range(Nz):

        # Assign nice names to primary variables
        h = sim.soln.h[:,:,ii]
        u = sim.soln.u[:,:,ii]
        v = sim.soln.v[:,:,ii]

        # Compute secondary varibles
        U = sim.avx_h(h)*u                  
        V = sim.avy_h(h)*v
        B = sim.gs[ii]*h + 0.5*(sim.avx_u(u**2) + sim.avy_v(v**2))
        q = (sim.ddx_v(v,dx) - sim.ddy_u(u,dy)  + sim.F)/(sim.avy_u(sim.avx_h(h)))  

        # Flux
        #sim.curr_flux.u[:,:,ii] =   sim.avy_v(q*sim.avx_v(V)) - sim.ddx_h(B,dx)
        #sim.curr_flux.v[:,:,ii] = - sim.avx_u(q*sim.avy_u(U)) - sim.ddy_h(B,dy)
        sim.curr_flux.u[:,:,ii] =   sim.avy_v(q)*sim.avy_v(sim.avx_v(V)) - sim.ddx_h(B,dx)
        sim.curr_flux.v[:,:,ii] = - sim.avx_u(q)*sim.avx_u(sim.avy_u(U)) - sim.ddy_h(B,dy)
        sim.curr_flux.h[:,:,ii] = - sim.ddx_u(U,dx) - sim.ddy_v(V,dy)
    return

def sadourny_sw_linear_flux(sim):

    Nx, Ny, Nz = sim.Nx, sim.Ny, sim.Nz
    dx, dy     = sim.dx[0], sim.dx[1]
    #ddx, ddy   = sim.ddx, sim.ddy
    #avx, avy   = sim.avx, sim.avy
    Hs         = sim.Hs[0]
    
    # Loop through each layer and compute the flux
    for ii in range(sim.Nz):

        # Assign nice names to primary variables
        h = sim.soln.h[:,:,ii]
        u = sim.soln.u[:,:,ii]
        v = sim.soln.v[:,:,ii]

        # Compute secondary varibles
        U = Hs*u             
        V = Hs*v
        q = sim.F/Hs
        B = sim.gs[ii]*h

        # Flux
        #sim.curr_flux.u[:,:,ii] =   sim.avy_v(q*sim.avx_v(V)) - sim.ddx_h(B,dx)
        #sim.curr_flux.v[:,:,ii] = - sim.avx_u(q*sim.avy_u(U)) - sim.ddy_h(B,dy)
        sim.curr_flux.u[:,:,ii] =   sim.avy_v(q)*sim.avy_v(sim.avx_v(V)) - sim.ddx_h(B,dx)
        sim.curr_flux.v[:,:,ii] = - sim.avx_u(q)*sim.avx_u(sim.avy_u(U)) - sim.ddy_h(B,dy)
        sim.curr_flux.h[:,:,ii] = - sim.ddx_u(U,dx) - sim.ddy_v(V,dy)

    return

def sadourny_sw(sim):

    # FJP: work on BCs
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

    sim.x_derivs = Diff.SADOURNY_x
    sim.y_derivs = Diff.SADOURNY_y

    if sim.dynamics == 'Nonlinear':
        sim.flux_function = sadourny_sw_flux
    elif sim.dynamics == 'Linear':
        sim.flux_function = sadourny_sw_linear_flux
    else:
        print "dynamics must be from the list: Nonlinear, Linear"
        sys.exit()
            
