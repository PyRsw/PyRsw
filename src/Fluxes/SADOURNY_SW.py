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
    ddx, ddy   = sim.ddx, sim.ddy
    avx, avy   = sim.avx, sim.avy
    
    # Loop through each layer and compute the flux
    for ii in range(Nz):

        # Assign nice names to primary variables
        h = sim.soln.h[:,:,ii]
        u = sim.soln.u[:,:,ii]
        v = sim.soln.v[:,:,ii]

        # FJP: try in the x-direction
        # FJP: try in x-y plane
        # FJP: work on BCx
        # Compute secondary varibles
        U = avx(h)*u                  
        V = avy(h)*v
        B = sim.gs[ii]*h[0:Nx,0:Ny] + 0.5*(avx(u**2)[:,0:Ny] + avy(np.concatenate([v[:,-1:]**2,v**2],axis=1)))
        q = (ddx(v,dx) - ddy(u,dy)  + sim.F)/(avy(avx(h)))  

        # Flux
        tmp = q*avx(V)
        sim.curr_flux.u[:,   0:Ny,ii] =   avy(np.concatenate([tmp[:,-1:],tmp],axis=1)) - ddx(B,dx)
        sim.curr_flux.v[0:Nx, :,  ii] = - avx(q*avy(U)) - ddy(np.concatenate([B,B[:,0:1]],axis=1),dy)
        sim.curr_flux.h[0:Nx,0:Ny,ii] = - ddx(U,dx)     - ddy(np.concatenate([V[:,-1:],V],axis=1),dy)

    return

def sadourny_sw_linear_flux(sim):

    Nx, Ny, Nz = sim.Nx, sim.Ny, sim.Nz
    dx, dy     = sim.dx[0], sim.dx[1]
    ddx, ddy   = sim.ddx, sim.ddy
    avx, avy   = sim.avx, sim.avy
    Hs         = sim.Hs[0]
    
    # Loop through each layer and compute the flux
    for ii in range(sim.Nz):

        # Assign nice names to primary variables
        h = sim.soln.h[:,:,ii]
        u = sim.soln.u[:,:,ii]
        v = sim.soln.v[:,:,ii]

        # Compute secondary varibles
        #   [U,V] = h[u,v]                  Transport velocities
        #   B = g*h + 0.5*(u**2 + v**2)     Bernoulli function
        #   q = (v_x - u_y + f)/h           Potential Vorticity
        U = Hs*u             
        V = Hs*v
        q = sim.F/Hs
        B = sim.gs[ii]*h[0:Nx,0:Ny]
        #FJP: maybe extend avy term so it's easier later?

        # Evolution Eqns:
        #	u_t =   (q*V^x)^y - d_x B
        #	v_t = - (q*U^y)^x - d_y B
        #	h_t = - H*u_x     - H*v_y
        # Flux
        tmp = q*avx(V)
        sim.curr_flux.u[:,   0:Ny,ii] =   avy(np.concatenate([tmp[:,-1:],tmp],axis=1)) - ddx(B,dx)
        sim.curr_flux.v[0:Nx, :,  ii] = - avx(q*avy(U)) - ddy(np.concatenate([B,B[:,0:1]],axis=1),dy)
        sim.curr_flux.h[0:Nx,0:Ny,ii] = - ddx(U,dx)     - ddy(np.concatenate([V[:,-1:],V],axis=1),dy)

    return

def sadourny_sw(sim):

    #FJP: impose ghost cells?
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
            
