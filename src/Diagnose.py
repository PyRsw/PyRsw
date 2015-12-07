# This file contains the diagnostics code
import numpy as np
import matplotlib.pyplot as plt

eps = 1e-10

# Initialize Diagnostics
def initialize_diagnostics(sim):
    sim.diag_times = [0.]
    sim.KEs = []
    sim.PEs = []
    sim.ENs = []
    sim.Ms  = []
    sim.next_diag_time = (np.floor(sim.time/sim.diagt)+1)*sim.diagt

    # Compute the initial values
    area = sim.dx[0]*sim.dx[1]
    KE = 0
    PE = 0
    EN = 0
    M  = 0

    for ii in range(sim.Nz):
        h   = sim.soln.h[:,:,ii] - sim.soln.h[:,:,ii+1]

        if sim.method == 'Spectral':
            v   = sim.soln.u[:,:,ii].copy()
            u   = sim.soln.v[:,:,ii].copy()
        else:
            v   = sim.soln.u[:,:,ii]/(eps + h)
            u   = sim.soln.v[:,:,ii]/(eps + h)
        dH2 = sim.soln.h[:,:,ii]**2 - sim.soln.h[:,:,ii+1]**2

        KE += 0.5*sim.rho[ii]*np.sum(np.ravel(area*h*(u**2 + v**2)))

        PE += area*0.5*sim.rho[ii]*sim.gs[0,ii]*np.sum(np.ravel(dH2))

        M += area*np.sum(np.ravel(h))*sim.rho[ii]

    sim.KEs += [KE]
    sim.PEs += [PE]
    sim.ENs += [EN]
    sim.Ms  += [M]

# Total KE
def compute_KE(sim):
    KE = 0
    area = sim.dx[0]*sim.dx[1]
    for ii in range(sim.Nz):
        h   = sim.soln.h[:,:,ii] - sim.soln.h[:,:,ii+1]
        if sim.method == 'Spectral':
            v   = sim.soln.u[:,:,ii].copy()
            u   = sim.soln.v[:,:,ii].copy()
        else:
            v   = sim.soln.u[:,:,ii]/(eps + h)
            u   = sim.soln.v[:,:,ii]/(eps + h)
        KE += 0.5*sim.rho[ii]*np.sum(np.ravel(area*h*(u**2 + v**2))) 
    return KE

# Total PE
def compute_PE(sim):
    PE = 0
    area = sim.dx[0]*sim.dx[1]
    for ii in range(sim.Nz):
        dH2 = sim.soln.h[:,:,ii]**2 - sim.soln.h[:,:,ii+1]**2
        PE += area*0.5*sim.rho[ii]*sim.gs[0,ii]*np.sum(np.ravel(dH2))
    return PE

# Compute enstrophy
def compute_enstrophy(sim):
    EN = 0
    for ii in range(sim.Nz-1,-1,-1):
        # Compute the mid-depth of the next layer
        h   = sim.soln.h[:,:,ii] - sim.soln.h[:,:,ii+1]
        if sim.method == 'Spectral':
            v   = sim.soln.u[:,:,ii].copy
            u   = sim.soln.v[:,:,ii].copy()
        else:
            v   = sim.soln.u[:,:,ii]/(eps + h)
            u   = sim.soln.v[:,:,ii]/(eps + h)
        q2  = (sim.dxp(v,sim.dx) - sim.dyp(u,sim.dx) + sim.f0)**2/h
        EN += np.sum(np.ravel(q2))
    return EN

def compute_mass(sim):
    M = 0.
    Area = sim.dx[0]*sim.dx[1]
    for ii in range(sim.Nz):
        h = sim.soln.h[:,:,ii] - sim.soln.h[:,:,ii+1]
        M += Area*np.sum(np.ravel(h))*sim.rho[ii]
    return M
    
# Update diagnostics
def update(sim):

    area = sim.dx[0]*sim.dx[1]
    KE = 0
    PE = 0
    EN = 0
    M  = 0

    for ii in range(sim.Nz):
        h   = sim.soln.h[:,:,ii] - sim.soln.h[:,:,ii+1]
        if sim.method == 'Spectral':
            v   = sim.soln.u[:,:,ii].copy()
            u   = sim.soln.v[:,:,ii].copy()
        else:
            v   = sim.soln.u[:,:,ii]/(eps + h)
            u   = sim.soln.v[:,:,ii]/(eps + h)
        dH2 = sim.soln.h[:,:,ii]**2 - sim.soln.h[:,:,ii+1]**2

        KE += area*0.5*sim.rho[ii]*np.sum(np.ravel(h*(u**2 + v**2)))

        PE += area*0.5*sim.rho[ii]*sim.gs[0,ii]*np.sum(np.ravel(dH2))

        M += area*np.sum(np.ravel(h))*sim.rho[ii]

    sim.KEs += [KE]
    sim.PEs += [PE]
    sim.ENs += [EN]
    sim.Ms  += [M]

    sim.diag_times += [sim.time]
    sim.next_diag_time += sim.diagt

def plot(sim):
    KE = np.array(sim.KEs)
    PE = np.array(sim.PEs)
    EN = np.array(sim.ENs)
    M  = np.array(sim.Ms)
    T  = np.array(sim.diag_times)

    fig = plt.figure()
            
    plt.subplot(2,2,1)
    plt.plot(T,KE - KE[0])
    plt.title('Tot. KE Dev.')
    plt.tight_layout()
    plt.locator_params(nbins=5)

    plt.subplot(2,2,2)
    plt.plot(T,PE - PE[0])
    plt.title('Tot. PE Dev.')
    plt.tight_layout()
    plt.locator_params(nbins=5)

    plt.subplot(2,2,3)
    plt.plot(T,(abs(KE) + abs(PE))/(KE[0]+PE[0]) - 1)
    plt.title('Rel. Energy Dev.')
    plt.tight_layout()
    plt.locator_params(nbins=5)

    plt.subplot(2,2,4)
    plt.plot(T,M/M[0] - 1.0)
    plt.title('Rel. Mass Dev.')
    plt.tight_layout()
    plt.locator_params(nbins=5)

    fig.tight_layout()
    fig.savefig('Outputs/{0:s}/diagnostics.pdf'.format(sim.run_name))
    
    return fig

def save(sim):
    if sim.output:
        fname = 'Outputs/{0:s}/diagnostics'.format(sim.run_name)
        np.savez_compressed(fname, KE = np.array(sim.KEs), 
                                   PE = np.array(sim.PEs),
                                   EN = np.array(sim.ENs),
                                   M  = np.array(sim.Ms))
