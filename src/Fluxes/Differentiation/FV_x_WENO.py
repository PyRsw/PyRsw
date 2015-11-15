import numpy as np
import sys

# Use a WENO 3-5 interpolation scheme
# Project onto the i+1/2 point
def ap(f):
    gam_p_L,gam_p_C,gam_p_R = 1./10, 6./10, 3./10
    gam_m_L,gam_m_C,gam_m_R = 3./10, 6./10, 1./10

    tmp1 = np.roll(f,2,0)
    tmp2 = np.roll(f,1,0)
    tmp3 = f.copy()

    # First stencil
    betaL = (13./12)*((tmp1 - 2.*tmp2 + tmp3)**2) + ( 1./4 )*((   tmp1 - 4.*tmp2 + 3.*tmp3)**2)
    pL =  (2./6)*tmp1 - (7./6)*tmp2  + (11./6)*tmp3
    mL = -(1./6)*tmp1 + (5./6)*tmp2  + (1./3)*tmp3

    # Second stencil
    tmp1 = np.roll(f,-1,0)
    betaC = (13./12)*((tmp2 - 2.*tmp3 + tmp1)**2) + ( 1./4 )*((   tmp2              - tmp1)**2)
    pC = -(1./6)*tmp2 + (5./6)*tmp3 + (2./6)*tmp1
    mC =  (1./3)*tmp2 + (5./6)*tmp3 - (1./6)*tmp1

    # Third stencil
    tmp2 = np.roll(f,-2,0)
    betaR = (13./12)*((tmp3 - 2.*tmp1 + tmp2)**2) + ( 1./4 )*((3.*tmp3 - 4.*tmp1    + tmp2)**2)
    pR = ( 1./3)*tmp3 + (5./6)*tmp1 - (1./6)*tmp2
    mR = (11./6)*tmp3 - (7./6)*tmp1 + (1./3)*tmp2

    eps = 1.e-6

    # Compute the modified coefficients
    wpL = gam_p_L/((eps+betaL)**2)
    wpC = gam_p_C/((eps+betaC)**2)
    wpR = gam_p_R/((eps+betaR)**2)
    wmL = gam_m_L/((eps+betaL)**2)
    wmC = gam_m_C/((eps+betaC)**2)
    wmR = gam_m_R/((eps+betaR)**2)

    sp  = wpL + wpC + wpR
    sm  = wmL + wmC + wmR

    wpL = wpL/sp
    wpC = wpC/sp
    wpR = wpR/sp
    wmL = wmL/sm
    wmC = wmC/sm
    wmR = wmR/sm

    ap = wpL*pL + wpC*pC + wpR*pR
    am = wmL*mL + wmC*mC + wmR*mR

    return am, ap

# Compute the mean derivative
# Remain at point i
def dp(f,dx):
    # Compute for downwinding
    gamma1,gamma2,gamma3 = 1./6,4./6,1./6

    tmp1 = np.roll(f,2,0)
    tmp2 = np.roll(f,1,0)
    tmp3 = f.copy()

    # First stencil
    beta1 = (13./12)*((tmp1 - 2.*tmp2 + tmp3)**2) + ( 1./4 )*((   tmp1 - 4.*tmp2 + 3.*tmp3)**2)
    diff1 = 0.5*(tmp1 - 4*tmp2 + 3*tmp3)/dx[0]

    # Second stencil
    tmp1 = np.roll(f,-1,0)
    beta2 = (13./12)*((tmp2 - 2.*tmp3 + tmp1)**2) + ( 1./4 )*((   tmp2              - tmp1)**2)
    diff2 = -0.5*(tmp2 - tmp1)/dx[0]

    # Third stencil
    tmp2 = np.roll(f,-2,0)
    beta3 = (13./12)*((tmp3 - 2.*tmp1 + tmp2)**2) + ( 1./4 )*((3.*tmp3 - 4.*tmp1    + tmp2)**2)
    diff3 = -0.5*(3*tmp3 - 4*tmp1 + tmp2)/dx[0]

    eps = 1.e-6

    # Compute the modified coefficients
    w1 = gamma1/((eps+beta1)**2)
    w2 = gamma2/((eps+beta2)**2)
    w3 = gamma3/((eps+beta3)**2)

    s  = w1 + w2 + w3

    w1 = w1/s
    w2 = w2/s
    w3 = w3/s

    diff = w1*diff1 + w2*diff2 + w3*diff3

    return diff

def DoNothing(f):
    return f,f
def DoNothing2(f,dx):
    return 0.*f

def FV_x_WENO(sim):
    if sim.Nx == 1:
        sim.axu = DoNothing
        sim.axv = DoNothing
        sim.axh = DoNothing

        sim.Dxu = DoNothing2
        sim.Dxv = DoNothing2
        sim.Dxh = DoNothing2
    else:
        if sim.BCx == 'periodic':
            try:
                import WENOx_periodic # try importing fortran code
                sim.axu = WENOx_periodic.ax
                sim.axv = WENOx_periodic.ax
                sim.axh = WENOx_periodic.ax

                sim.Dxu = WENOx_periodic.ddx
                sim.Dxv = WENOx_periodic.ddx
                sim.Dxh = WENOx_periodic.ddx
                sim.Dx = WENOx_periodic.ddx
                print('Using fortran in x')
            except:
                print('Unable to load fortran for x interpolation with wall boundary conditions.')
                sys.exit()

                sim.axu = ap
                sim.axv = ap
                sim.axh = ap

                sim.Dxu = dp
                sim.Dxv = dp
                sim.Dxh = dp
                print('Not using fortran in x')
        if sim.BCx == 'wall':
            try:
                import WENOx_dirichlet # try importing fortran code
                import WENOx_neumann # try importing fortran code
                sim.axu = WENOx_dirichlet.ax
                sim.axv = WENOx_neumann.ax
                sim.axh = WENOx_neumann.ax

                sim.Dxu = WENOx_dirichlet.ddx
                sim.Dxv = WENOx_neumann.ddx
                sim.Dxh = WENOx_neumann.ddx
                sim.Dx = WENOx_neumann.ddx
                print('Using fortran in x')
            except:
                print('Unable to load fortran for x interpolation with wall boundary conditions.')
                sys.exit()
