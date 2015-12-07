import numpy as np
import sys

def RK4(sim):

    # y^{n+1} = y^n + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
    # k1 = flux(y^n)
    # k2 = flux(y^n + (dt/2)*k1)
    # k3 = flux(y^n + (dt/2)*k2)
    # k4 = flux(y^n + dt*k3)

    # Compute k1
    flux_u, flux_v, flux_h = sim.flux()
    src_u,  src_v,  src_h  = sim.source()
    k1_u = flux_u + src_u
    k1_v = flux_v + src_v
    k1_h = flux_h + src_h

    # Compute k2
    sim.soln.u += k1_u*sim.dt/2.
    sim.soln.v += k1_v*sim.dt/2.
    sim.soln.h += k1_h*sim.dt/2.
    flux_u, flux_v, flux_h = sim.flux()
    src_u,  src_v,  src_h  = sim.source()
    k2_u = flux_u + src_u
    k2_v = flux_v + src_v
    k2_h = flux_h + src_h

    # Compute k3
    sim.soln.u += -k1_u*sim.dt/2. + k2_u*sim.dt/2.
    sim.soln.v += -k1_v*sim.dt/2. + k2_v*sim.dt/2.
    sim.soln.h += -k1_h*sim.dt/2. + k2_h*sim.dt/2.
    flux_u, flux_v, flux_h = sim.flux()
    src_u,  src_v,  src_h  = sim.source()
    k3_u = flux_u + src_u
    k3_v = flux_v + src_v
    k3_h = flux_h + src_h

    # Compute k4
    sim.soln.u += -k2_u*sim.dt/2. + k3_u*sim.dt
    sim.soln.v += -k2_v*sim.dt/2. + k3_v*sim.dt
    sim.soln.h += -k2_h*sim.dt/2. + k3_h*sim.dt
    flux_u, flux_v, flux_h = sim.flux()
    src_u,  src_v,  src_h  = sim.source()
    k4_u = flux_u + src_u
    k4_v = flux_v + src_v
    k4_h = flux_h + src_h

    sim.soln.u += -k3_u*sim.dt + (sim.dt/6.)*(k1_u + 2.*k2_u + 2.*k3_u + k4_u)
    sim.soln.v += -k3_v*sim.dt + (sim.dt/6.)*(k1_v + 2.*k2_v + 2.*k3_v + k4_v)
    sim.soln.h += -k3_h*sim.dt + (sim.dt/6.)*(k1_h + 2.*k2_h + 2.*k3_h + k4_h)
