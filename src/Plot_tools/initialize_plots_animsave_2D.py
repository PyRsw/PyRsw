# Initialize plot objects for anim or save
# Assume that hte field is 2-dimensional

import numpy as np
import matplotlib.pyplot as plt

from update_anim_2D import update_anim_2D
from update_save_2D import update_save_2D

def initialize_plots_animsave_2D(sim):

    figs = []
    Qs   = []
    ttls = []

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
                to_plot = sim.soln.u[0:sim.Nx,0:sim.Ny,L]
            elif var == 'v':
                ttl = fig.suptitle('Meridional Velocity : t = 0')
                to_plot = sim.soln.v[0:sim.Nx,0:sim.Ny,L]
            elif var == 'h':
                ttl = fig.suptitle('Free Surface Displacement : t = 0')
                to_plot = sim.soln.h[0:sim.Nx,0:sim.Ny,L] - sim.Hs[L]
            elif var == 'vort':
                to_plot =     sim.ddx_v(sim.soln.v[0:sim.Nx,0:sim.Ny,L],sim) \
                            - sim.ddy_u(sim.soln.u[0:sim.Nx,0:sim.Ny,L],sim)
                if sim.f0 != 0:
                    ttl = fig.suptitle('Vorticity / f_0 : t = 0')
                    to_plot *= 1./sim.f0
                else:   
                    ttl = fig.suptitle('Vorticity : t = 0')
            elif var == 'div':
                h = sim.soln.h[:,:,L] 
                to_plot =     sim.ddx_u(h*sim.soln.u[0:sim.Nx,0:sim.Ny,L],sim) \
                            + sim.ddy_v(h*sim.soln.v[0:sim.Nx,0:sim.Ny,L],sim)
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

            Q = plt.pcolormesh(sim.X/1e3, sim.Y/1e3, to_plot, cmap=sim.cmap, 
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
        plt.ioff()
        plt.pause(0.01)

    if sim.animate == 'Anim':
        plt.ion()
        plt.pause(0.01)
        plt.draw()

    sim.figs = figs
    sim.Qs = Qs
    sim.ttls = ttls

