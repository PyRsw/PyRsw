import numpy as np
from Euler import Euler
import matplotlib.pyplot as plt

def AB2(sim):
    if sim.nfluxes < 1:
        sim.nfluxes = 1

    if len(sim.fluxes.u) == 0:

        Euler(sim)

    elif len(sim.fluxes.u) == 1:
        
        sim.flux()

        w1 = sim.dt*(1. + 0.5*sim.dt/sim.dts[0])
        w2 = -0.5*sim.dt**2/sim.dts[0]
        
        sim.soln.u += w1*sim.curr_flux.u + w2*sim.fluxes.u[0]
        sim.soln.v += w1*sim.curr_flux.v + w2*sim.fluxes.v[0]
        sim.soln.h += w1*sim.curr_flux.h + w2*sim.fluxes.h[0]
        
        if sim.nfluxes == 1:
            sim.fluxes.u = [sim.curr_flux.u]
            sim.fluxes.v = [sim.curr_flux.v]
            sim.fluxes.h = [sim.curr_flux.h]
            sim.dts    = [sim.dt]
        else:
            sim.fluxes.u = [sim.curr_flux.u] + sim.fluxes.u
            sim.fluxes.v = [sim.curr_flux.v] + sim.fluxes.v
            sim.fluxes.h = [sim.curr_flux.h] + sim.fluxes.h
            sim.dts = [sim.dt] + sim.dts
