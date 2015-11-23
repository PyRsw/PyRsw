import numpy as np

def Euler(sim):

    # Compute the flux
    sim.flux()
    
    # Evolve the system
    sim.soln.u += sim.curr_flux.u*sim.dt
    sim.soln.v += sim.curr_flux.v*sim.dt
    sim.soln.h += sim.curr_flux.h*sim.dt
    
    # Store the appropriate history
    if sim.nfluxes > 0:
        sim.fluxes.u = [sim.curr_flux.u.copy()]
        sim.fluxes.v = [sim.curr_flux.v.copy()]
        sim.fluxes.h = [sim.curr_flux.h.copy()]
        sim.dts    = [sim.dt]
