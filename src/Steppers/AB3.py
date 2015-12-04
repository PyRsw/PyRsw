import numpy as np
from Euler import Euler
from AB2 import AB2
import matplotlib.pyplot as plt


def AB3(sim):
    if sim.nfluxes < 2:
        sim.nfluxes = 2

    if len(sim.fluxes.u) == 0:

        Euler(sim)

    elif len(sim.fluxes.u) == 1:

        AB2(sim)

    elif len(sim.fluxes.u) == 2:

        # Compute the fluxes
        sim.flux()

        # Compute the weights for the adaptive time-stepped
        # Adams-Bashforth 3 scheme
        tp = sim.time + sim.dt
        tn = sim.time

        a  = sim.time - sim.dts[0]
        b  = sim.time - sim.dts[0] - sim.dts[1]
        ts = sim.time

        w0 = (   (1./3)*(tp**3 - tn**3) \
                -   0.5*(tp**2 - tn**2)*(a+b) \
                +   a*b*(tp    - tn))   \
             /((ts-a)*(ts-b))

        a  = sim.time
        b  = sim.time - sim.dts[0] - sim.dts[1]
        ts = sim.time - sim.dts[0]

        w1 = (   (1./3)*(tp**3 - tn**3) \
                -   0.5*(tp**2 - tn**2)*(a+b) \
                +   a*b*(tp    - tn))   \
             /((ts-a)*(ts-b))

        a  = sim.time
        b  = sim.time - sim.dts[0]
        ts = sim.time - sim.dts[0] - sim.dts[1]

        w2 = (   (1./3)*(tp**3 - tn**3) \
                -   0.5*(tp**2 - tn**2)*(a+b) \
                +   a*b*(tp    - tn))   \
             /((ts-a)*(ts-b))

        # Evolve the system
        sim.soln.u += w0*sim.curr_flux.u + w1*sim.fluxes.u[0] + w2*sim.fluxes.u[1]
        sim.soln.v += w0*sim.curr_flux.v + w1*sim.fluxes.v[0] + w2*sim.fluxes.v[1]
        sim.soln.h += w0*sim.curr_flux.h + w1*sim.fluxes.h[0] + w2*sim.fluxes.h[1]
        
        # Store the appropriate histories.
        if sim.nfluxes == 2:
            sim.fluxes.u = [sim.curr_flux.u.copy(), sim.fluxes.u[0].copy()]
            sim.fluxes.v = [sim.curr_flux.v.copy(), sim.fluxes.v[0].copy()]
            sim.fluxes.h = [sim.curr_flux.h.copy(), sim.fluxes.h[0].copy()]
            sim.dts    = [sim.dt, sim.dts[0]]
        else:
            sim.fluxes.u = [sim.curr_flux.u.copy()] + sim.fluxes.u
            sim.fluxes.v = [sim.curr_flux.v.copy()] + sim.fluxes.v
            sim.fluxes.h = [sim.curr_flux.h.copy()] + sim.fluxes.h
            sim.dts = [sim.dt] + sim.dts
