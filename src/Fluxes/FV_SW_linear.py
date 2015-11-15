import numpy as np
import matplotlib.pyplot as plt
import sys

# Compute the flux terms
def fv_sw_linear_flux(sim):

    # Some storage fields
    the_flux = np.zeros(sim.sol.shape)
    hs = np.zeros((sim.Nx,sim.Ny,sim.Nz+1))

    # Difference in layer deformations gives layer depths
    hs[:,:,:sim.Nz] = sim.sol[sim.Ih,:,:,:sim.Nz] - sim.sol[sim.Ih,:,:,1:]

    # Loop through each layer and compute flux
    for ii in range(sim.Nz):
       
        h  = hs[:,:,ii]
        uh = sim.sol[sim.Iu,:,:,ii]
        vh = sim.sol[sim.Iv,:,:,ii]

        # Compute the interpolations
        hmx,hpx = sim.ax(h)
        hmy,hpy = sim.ay(h)

        umx,upx = sim.ax(uh)
        umy,upy = sim.ay(uh)
#        umx = umx/hmx; upx = upx/hpx
#        umy = umy/hmy; upy = upy/hpy

        vmx,vpx = sim.ax(vh)
        vmy,vpy = sim.ay(vh)
#        vmx = vmx/hmx; vpx = vpx/hpx
#        vmy = vmy/hmy; vpy = vpy/hpy

        # Compute the winds
        EastWest = uh.copy()
        NorthSouth = vh.copy()

        EW = EastWest >= 0. # East Wind
        WW = EastWest < 0. # West Wind
        NW = NorthSouth >= 0. # North Wind
        SW = NorthSouth < 0. # South Wind

        # Compute fluxes

        ##
        ## u flux
        ## (0.5*g*h*h)_x 
        ##

        #FUp = 0.5*sim.g*hpx*hpx
        FUp = sim.g*hpx
        FUm = np.roll(FUp,1,0)
        #FDp = 0.5*sim.g*hmx*hmx
        FDp = sim.g*hmx
        FDm = np.roll(FDp,-1,0)

        Fxp = WW*FUp - EW*FDp 
        Fxm = WW*FUm - EW*FDm 

        the_flux[sim.Iu,:,:,ii] =  -(Fxp - Fxm)/sim.dx[0] 
        
        ##
        ## v flux
        ## (0.5*g*h*h)_y
        ##

        #FUp = 0.5*sim.g*hpy*hpy
        FUp = sim.g*hpy
        FUm = np.roll(FUp,1,1)
        #FDp = 0.5*sim.g*hmy*hmy
        FDp = sim.g*hmy
        FDm = np.roll(FDp,-1,1)

        Fyp = SW*FUp - NW*FDp 
        Fym = SW*FUm - NW*FDm 

        the_flux[sim.Iv,:,:,ii] =  - (Fyp - Fym)/sim.dx[1]

        ##
        ## h flux
        ## H0(u)_x + H0(v)_y
        ##

        # H0(u)_x 

        FUp = sim.Hs[ii]*upx
        FUm = np.roll(FUp,1,0)
        FDp = sim.Hs[ii]*umx
        FDm = np.roll(FDp,-1,0)

        Fxp = EW*FUp - WW*FDp
        Fxm = EW*FUm - WW*FDm

        # H0(v)_y 

        FUp = sim.Hs[ii]*vpy
        FUm = np.roll(FUp,1,1)
        FDp = sim.Hs[ii]*vmy
        FDm = np.roll(FDp,-1,1)

        Fyp = NW*FUp - SW*FDp
        Fym = NW*FUm - SW*FDm

        the_flux[sim.Ih,:,:,ii] =  -(Fxp - Fxm)/sim.dx[0] - (Fyp - Fym)/sim.dx[1]

    return the_flux

# Deal with the source terms
def fv_sw_linear_source(sim):
    
    # Some storage fields
    the_source = np.zeros(sim.sol.shape)
    h = np.zeros((sim.Nx,sim.Ny,sim.Nz+1))
    p = np.zeros((sim.Nx,sim.Ny,sim.Nz))

    # Difference in layer deformations gives layer depths
    h[:,:,:sim.Nz] = sim.sol[sim.Ih,:,:,:sim.Nz] - sim.sol[sim.Ih,:,:,1:]

    # Loop through each layer and compute flux

    # For layer i, the source terms look like:
    #  -g * h_i / rho_i * sum(rho_j * (d/dx)h_j) - g * h_j * (d/dx)*eta_{i+1}
    # Since the sum is built incrementally, I'll store one variable and build it as we go
    Sx = np.zeros(sim.sol.shape)
    Sy = np.zeros(sim.sol.shape)

    # Deal with first layer
    the_source[sim.Iu,:,:,0] = -sim.g*h[:,:,0]*sim.Dx(sim.sol[sim.Ih,:,:,1],sim.dx) \
                                    + sim.f0*sim.sol[sim.Iv,:,:,0]
    the_source[sim.Iv,:,:,0] = -sim.g*h[:,:,0]*sim.Dy(sim.sol[sim.Ih,:,:,1],sim.dx) \
                                    - sim.f0*sim.sol[sim.Iu,:,:,0]

    # Deal with rest of layers
    for ii in range(1,sim.Nz):
        Sx += sim.rho[ii-1]*sim.Dx(h[:,:,ii-1],sim.dx)
        Sy += sim.rho[ii-1]*sim.Dy(h[:,:,ii-1],sim.dx)

        # Deal with sum term
        the_source[sim.Iu,:,:,ii] += -(sim.g/sim.rho[ii])*h[:,:,ii]*Sx
        the_source[sim.Iv,:,:,ii] += -(sim.g/sim.rho[ii])*h[:,:,ii]*Sy

        # Deal with Coriolis terms
        the_source[sim.Iu,:,:,ii] +=  sim.f0*sim.sol[sim.Iv,:,:,ii]
        the_source[sim.Iv,:,:,ii] += -sim.f0*sim.sol[sim.Iu,:,:,ii]

        # Deal with final term
        the_source[sim.Iu,:,:,ii] += -sim.g*h[:,:,ii]*sim.Dx(sim.sol[sim.Ih,:,:,ii+1],sim.dx)
        the_source[sim.Iv,:,:,ii] += -sim.g*h[:,:,ii]*sim.Dy(sim.sol[sim.Ih,:,:,ii+1],sim.dx)

    return 0.*the_source

