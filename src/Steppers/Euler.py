import numpy as np

def Euler(sim):
    sim.flux()
    
    sim.soln.u += sim.curr_flux.u*sim.dt
    sim.soln.v += sim.curr_flux.v*sim.dt
    sim.soln.h += sim.curr_flux.h*sim.dt
    
    if sim.nfluxes > 0:
        sim.fluxes.u = [sim.curr_flux.u]
        sim.fluxes.v = [sim.curr_flux.v]
        sim.fluxes.h = [sim.curr_flux.h]
        sim.dts    = [sim.dt]
