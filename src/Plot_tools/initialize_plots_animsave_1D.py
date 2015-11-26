# Initialize plot objects for anim or save
# Assume that the field is 1-dimensional

import matplotlib.pyplot as plt
import matplotlib.animation as anim
import numpy as np

from update_anim_1D import update_anim_1D
from update_save_1D import update_save_1D

def initialize_plots_animsave_1D(sim):

    fig = plt.figure()
    sim.fig = fig
    fig.suptitle('t = 0')

    if sim.Nx > 1:
        x = sim.x/1e3
        xlabel = 'x (km)'
    else:
        x = sim.y/1e3
        xlabel = 'y (km)'

    Qs  = []
    axs = []

    # Loop through each element of plot_vars
    # Each field with be given its own subplot
    for var_cnt in range(len(sim.plot_vars)):

        var = sim.plot_vars[var_cnt]
        Qs  += [[]]
        axs += [[]]
    

        # Plot the data
        plt.subplot(len(sim.plot_vars),1,var_cnt+1)
        for L in range(sim.Nz):

            if var == 'u':
                to_plot = sim.soln.u[:,:,L].ravel()
            elif var == 'v':
                to_plot = sim.soln.v[:,:,L].ravel()
            elif var == 'h':
                to_plot = sim.soln.h[:,:,L].ravel() - sim.Hs[L]
            elif var == 'vort':
                to_plot = sim.ddx_v(sim.soln.v[:,:,L],sim) \
                        - sim.ddy_u(sim.soln.u[:,:,L],sim)
                to_plot = to_plot.ravel()
                if sim.f0 != 0:
                    to_plot *= 1./sim.f0
            elif var == 'div':
                h = sim.soln.h[:,:,L].ravel() 
                to_plot = sim.ddx_u(h*sim.soln.u[:,:,L],sim) \
                        + sim.ddy_v(h*sim.soln.v[:,:,L],sim)
                if sim.f0 != 0:
                    to_plot *= 1./sim.f0
                to_plot = to_plot.ravel()
            
            l, = plt.plot(x, to_plot, linewidth=2)

            # Has the user specified plot limits?
            if len(sim.ylims[var_cnt]) == 2:
                plt.ylim(sim.ylims[var_cnt])
            else:
                tmp = plt.gca().get_ylim()
                plt.ylim([-np.max(np.abs(tmp)), np.max(np.abs(tmp))]);

            plt.ylabel(var)

            # Store the plot for update purposes
            axs[var_cnt] += [plt.gca()]
            Qs[var_cnt] += [l]

        plt.xlabel(xlabel)


    if sim.animate == 'Anim':
        sim.update_plots = update_anim_1D
    elif sim.animate == 'Save':
        sim.update_plots = update_save_1D
        plt.pause(0.01)

    if sim.animate == 'Anim':
        plt.ion()
        plt.pause(0.01)
        plt.draw()
        
    sim.Qs  = Qs
    sim.axs = axs

