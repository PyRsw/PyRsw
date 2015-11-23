# Contains plotting commands
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import numpy as np

# Helper for formatting time strings
def smart_time(t):
    tstr = 't = '

    if t < 2*60.:
        tstr += '{0:.4f} sec'.format(t)
    elif t < 90.*60:
        tstr += '{0:.4f} min'.format(t/60)
    elif t < 48.*60*60:
        tstr += '{0:.4f} hrs'.format(t/(60*60))
    else:
        tstr += '{0:.4f} day'.format(t/(24*60*60))

    return tstr

# Update plot objects if animating
def update_anim(sim):

    sim.fig.suptitle(smart_time(sim.time))

    if sim.Nx > 1 and sim.Ny > 1:
        for L in range(sim.Nz):
            sim.Qs[L].set_array(np.ravel(sim.soln.u[:sim.Nx-1,:sim.Ny-1,L]))
            sim.Qs[L].changed()
    else:
        # Update u
        for L in range(sim.Nz):
            sim.Qs[0][L].set_ydata(sim.soln.u[:,:,L])
        sim.axs[0].relim()
        sim.axs[0].autoscale_view()

        # Update v
        for L in range(sim.Nz):
            sim.Qs[1][L].set_ydata(sim.soln.v[:,:,L])
        sim.axs[1].relim()
        sim.axs[1].autoscale_view()

        # Update h
        for L in range(sim.Nz):
            sim.Qs[2][L].set_ydata(sim.soln.h[:,:,L] - np.sum(sim.Hs[L:]))
        sim.axs[2].relim()
        sim.axs[2].autoscale_view()

    plt.pause(0.01) 
    plt.draw()

# Create Hovmuller plot at end of anim simulation
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


# Update plot objects if saving
def update_save(sim):

    sim.fig.suptitle(smart_time(sim.time))

    if sim.Nx > 1 and sim.Ny > 1:
        for L in range(sim.Nz):
            sim.Qs[L].set_array(np.ravel(sim.soln.u[:sim.Nx-1,:sim.Ny-1,L].T))
            sim.Qs[L].changed()
    else:
        # Update u
        for L in range(sim.Nz):
            sim.Qs[0][L].set_ydata(sim.soln.u[:,:,L])
        sim.axs[0].relim()
        sim.axs[0].autoscale_view()

        # Update v
        for L in range(sim.Nz):
            sim.Qs[1][L].set_ydata(sim.soln.v[:,:,L])
        sim.axs[1].relim()
        sim.axs[1].autoscale_view()

        # Update h
        for L in range(sim.Nz):
            sim.Qs[2][L].set_ydata(sim.soln.h[:,:,L] - np.sum(sim.Hs[L:]))
        sim.axs[2].relim()
        sim.axs[2].autoscale_view()

    plt.draw()

    sim.fig.savefig('Frames/{0:05d}.png'.format(sim.frame_count))
    sim.frame_count += 1


def update_hov(sim):

    sim.fig.suptitle(smart_time(sim.time))
    if sim.Nx > 1:
        sim.hov_h[:,:,sim.hov_count] = sim.soln.h[:,0,:-1]
        x = sim.x/1e3
    else:
        sim.hov_h[:,:,sim.hov_count] = sim.soln.h[0,:,:-1] 
        x = sim.y/1e3
    sim.hov_count += 1
    

# Initialize plot objects for anim or save
def initialize_plots_animsave(sim):

    fig = plt.figure()
    sim.fig = fig
    fig.suptitle('t = 0')

    if sim.Nx > 1 and sim.Ny > 1:
        x = sim.X/1e3
        y = sim.Y/1e3
        Qs  = []
        axs = []
        for L in range(sim.Nz):
            plt.subplot(1,sim.Nz,L+1)
            axs += [plt.gca()]
            Qs += [plt.pcolormesh(x,y,sim.soln.u[:,:,L], cmap='viridis')]
            plt.colorbar()
            try:
                plt.contour(x,y,sim.soln.u[:,:,-1])
            except:
                pass
            plt.axis('tight')
            plt.gca().set_aspect('equal')
            plt.pause(2.0)
            
    else:
        Qs  = [[],[],[]]
        axs = []

        # Plot u
        plt.subplot(3,1,1)
        axs += [plt.gca()]
        for L in range(sim.Nz):
            if sim.Nx > 1:
                x = sim.x/1e3
            if sim.Ny > 1:
                x = sim.y/1e3
            l, = plt.plot(x,sim.soln.u[:,:,L].ravel(), linewidth=2)
            if len(sim.ylims[0]) == 2:
                plt.ylim(sim.ylims[0])
            plt.ylabel('u')
            plt.ylim([-0.2, 0.2])
            Qs[0] += [l]

        # Plot v
        plt.subplot(3,1,2)
        axs += [plt.gca()]
        for L in range(sim.Nz):
            if sim.Nx > 1:
                x = sim.x/1e3
            if sim.Ny > 1:
                x = sim.y/1e3
            l, = plt.plot(x,sim.soln.v[:,:,L].ravel(), linewidth=2)
            if len(sim.ylims[1]) == 2:
                plt.ylim(sim.ylims[1])
            plt.ylabel('v')
            plt.ylim([-0.2, 0.2]);
            Qs[1] += [l]

        # Plot h
        plt.subplot(3,1,3)
        axs += [plt.gca()]
        for L in range(sim.Nz):
            if sim.Nx > 1:
                x = sim.x/1e3
            if sim.Ny > 1:
                x = sim.y/1e3
            l, = plt.plot(x,sim.soln.h[:,:,L].ravel() - np.sum(sim.Hs[L:]), linewidth=2)
            if len(sim.ylims[2]) == 2:
                plt.ylim(sim.ylims[2])
            plt.ylim([-0.2, 1.2]);
            plt.ylabel('eta')
            Qs[2] += [l]

        plt.plot(x,sim.soln.h[:,:,-1].ravel(),'k')

    if sim.animate == 'Anim':
        sim.update_plots = update_anim
    elif sim.animate == 'Save':
        sim.update_plots = update_save

    if sim.animate == 'Anim':
        plt.ion()
        plt.pause(0.01)
        plt.draw()
        
    sim.Qs  = Qs
    sim.axs = axs


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
        print(x.shape,t.shape,sim.hov_h.shape)
        plt.pcolormesh(x,t,-1000*np.ones(sim.hov_h[:,L,:].T.shape),vmin = 0, vmax = 1, cmap = 'cubehelix')
        #Q = plt.pcolormesh(x,t,sim.hov_h[:,L,:].T - np.sum(sim.Hs[L:]),cmap=sim.cmap)
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


# Finalize
def end_movie(sim):
    pass
