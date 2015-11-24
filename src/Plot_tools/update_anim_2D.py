# Update plot objects if animating
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import numpy as np
from smart_time import smart_time

def update_anim_2D(sim):

    for var_cnt in range(len(sim.plot_vars)):

        var = sim.plot_vars[var_cnt]
        sim.figs[var_cnt].suptitle('{0:s} : {1:s}'.format(var,smart_time(sim.time)))

        for L in range(sim.Nz):

            if var == 'u':
                to_plot = sim.soln.u[:,:,L]
            elif var == 'v':
                to_plot = sim.soln.v[:,:,L]
            elif var == 'h':
                to_plot = sim.soln.h[:,:,L] - sim.Hs[L]
            elif var == 'vort':
                to_plot =     sim.ddx_v(sim.soln.v[:,:,L],sim) \
                            - sim.ddy_u(sim.soln.u[:,:,L],sim)

            sim.Qs[var_cnt][L].set_array(to_plot.ravel())
            sim.Qs[var_cnt][L].changed()

        plt.pause(0.01) 
    plt.draw()

