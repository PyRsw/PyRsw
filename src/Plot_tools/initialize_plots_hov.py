import matplotlib.pyplot as plt
import numpy as np

# Initialize plot objects for hov
def initialize_plots_hov(sim):

    fig = plt.figure()
    sim.fig = fig
    fig.suptitle('t = 0')

    Qs = []

    # Plot h
    t = np.arange(0,sim.end_time+sim.plott,sim.plott)/86400.

    for L in range(sim.Nz):
        plt.subplot(sim.Nz,1,L+1)
        if sim.Nx > 1:
            x = sim.x/1e3
            plt.xlabel('x (km)')
        if sim.Ny > 1:
            x = sim.y/1e3
            plt.xlabel('y (km)')
        Q = plt.pcolormesh(x,t,sim.hov_h[:,L,:].T - np.sum(sim.Hs[L:]),cmap=sim.cmap)
        sim.vmin = range(sim.Nz)
        sim.vmax = range(sim.Nz)
        if len(sim.ylims[2]) == 2:
            sim.vmin[L] = sim.ylims[2][0]
            sim.vmax[L] = sim.ylims[2][1]
        else:
            tmp = np.max(np.abs(sim.hov_h[:,L,0] - np.sum(sim.Hs[L:])))
            sim.vmin[L] = -tmp
            sim.vmax[L] =  tmp
        Q = plt.pcolormesh(x,t[:2],sim.hov_h[:,L,:2].T - np.sum(sim.Hs[L:]),cmap=sim.cmap, 
                vmin = sim.vmin[L], vmax = sim.vmax[L])
        plt.ylim((t[0],t[-1]))
        plt.axis('tight')
        plt.title('eta')
        plt.colorbar()
        plt.ylabel('Time (days)')
        Qs += [Q]

    sim.update_plots = update_hov

    sim.Qs  = Qs

    plt.ion()
    plt.pause(0.01)
    plt.draw()

