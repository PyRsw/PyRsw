# Create Hovmuller plot at end of anim simulation
import matplotlib.pyplot as plt
import numpy as np

def plot_hov(sim):
    plt.figure()
    t = np.arange(0,sim.end_time+sim.plott,sim.plott)/86400.

    if sim.Ny==1:
        x = sim.x/1e3
    elif sim.Nx == 1:
        x = sim.y/1e3

    for L in range(sim.Nz):
        field = sim.hov_h[:,0,:].T - np.sum(sim.Hs[L:])
        cv = np.max(np.abs(field.ravel()))
        plt.subplot(sim.Nz,1,L+1)
        plt.pcolormesh(x,t, field,
            cmap=sim.cmap, vmin = -cv, vmax = cv)
        plt.axis('tight')
        plt.title(r"$\mathrm{Hovm{\"o}ller} \; \mathrm{Plot} \; \mathrm{of} \; \eta$", fontsize = 16)
        if sim.Nx > 1:
            plt.xlabel(r"$\mathrm{x} \; \mathrm{(km)}$", fontsize=14)
        else:
            plt.xlabel(r"$\mathrm{y} \; \mathrm{(km)}$", fontsize=14)
        plt.ylabel(r"$\mathrm{Time} \; \mathrm{(days)}$", fontsize=14)
        plt.colorbar()

