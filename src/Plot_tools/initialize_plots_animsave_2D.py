# Initialize plot objects for anim or save
# Assume that hte field is 2-dimensional

import matplotlib.pyplot as plt
import matplotlib.animation as anim
import numpy as np

from update_anim_2D import update_anim_2D
from update_save_2D import update_save_2D

def initialize_plots_animsave_2D(sim):

    figs = []
    Qs   = []
    ttls = []

    # Create the appropriate grids:
    # extend by a point to handle the way
    # that pcolor removes boundaries
    x = np.zeros((sim.Nx+1))
    x[1:] = sim.x + (sim.x[1] - sim.x[0])/2.
    x[0]  = sim.x[0] - (sim.x[1] - sim.x[0])/2.

    y = np.zeros((sim.Ny+1))
    y[1:] = sim.y + (sim.y[1] - sim.y[0])/2.
    y[0]  = sim.y[0] - (sim.y[1] - sim.y[0])/2.

    X,Y = np.meshgrid(x,y,indexing='ij')

    # Loop through each element of plot_vars
    # Each field will be given its own figure
    for var_cnt in range(len(sim.plot_vars)):

        Qs   += [[]]
        ttls += [[]]

        var = sim.plot_vars[var_cnt]

        fig = plt.figure()
        figs += [fig]

        # Plot the data
        for L in range(sim.Nz):

            plt.subplot(1,sim.Nz,L+1)

            if var == 'u':
                ttl = fig.suptitle('Zonal Velocity : t = 0')
                to_plot = sim.soln.u[:,:,L]
            elif var == 'v':
                ttl = fig.suptitle('Meridional Velocity : t = 0')
                to_plot = sim.soln.v[:,:,L]
            elif var == 'h':
                ttl = fig.suptitle('Free Surface Displacement : t = 0')
                to_plot = sim.soln.h[:,:,L] - sim.Hs[L]
            elif var == 'vort':
                to_plot =     sim.ddx_v(sim.soln.v[:,:,L],sim) \
                            - sim.ddy_u(sim.soln.u[:,:,L],sim)
                if sim.f0 != 0:
                    ttl = fig.suptitle('Vorticity / f_0 : t = 0')
                    to_plot *= 1./sim.f0
                else:   
                    ttl = fig.suptitle('Vorticity : t = 0')
            elif var == 'div':
                h = sim.soln.u[:,:,L] + sim.Hs[L]
                to_plot =     sim.ddx_u(h*sim.soln.u[:,:,L],sim) \
                            + sim.ddy_v(h*sim.soln.v[:,:,L],sim)
                if sim.f0 != 0:
                    ttl = fig.suptitle('Divergence of mass-flux / f_0 : t = 0')
                    to_plot *= 1./sim.f0
                else:   
                    ttl = fig.suptitle('Divergence of mass-flux : t = 0')

            # Has the user specified plot limits?
            if len(sim.clims[var_cnt]) == 2:
                vmin = sim.clims[var_cnt][0]
                vmax = sim.clims[var_cnt][1]
            else:
                cv = np.max(np.abs(to_plot.ravel()))
                vmin = -cv
                vmax =  cv

            Q = plt.pcolormesh(X/1e3, Y/1e3, to_plot, cmap=sim.cmap, 
                        vmin = vmin, vmax = vmax)
            Qs[var_cnt] += [Q]
            ttls[var_cnt] += [ttl]

            plt.colorbar()

            plt.axis('tight')
            plt.gca().set_aspect('equal')
            

    if sim.animate == 'Anim':
        sim.update_plots = update_anim_2D
    elif sim.animate == 'Save':
        sim.update_plots = update_save_2D
        plt.pause(0.01)

    if sim.animate == 'Anim':
        plt.ion()
        plt.pause(0.01)
        plt.draw()
        
    sim.figs = figs
    sim.Qs = Qs
    sim.ttls = ttls

