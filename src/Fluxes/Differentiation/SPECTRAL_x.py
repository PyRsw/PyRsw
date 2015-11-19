import numpy as np
import sys
 
try:
    from pyfftw.interfaces.scipy_fftpack import fft, ifft, fftfreq, fftn, ifftn
except:
    from scipy.fftpack import fft, ifft, fftfreq, fftn, ifftn

def ddx_period(f,ik):

    df = np.real(ifft(ik*fft(f,axis=0),axis=0))

    return df

def ddx_even(f,ik):

    N = f.shape[0]
    fe = np.concatenate([f,f[::-1,:]],axis=0)
    df = np.real(ifft(ik*fft(fe,axis=0),axis=0))
    df = df[:N,:]
    
    return df

def ddx_odd(f,ik):

    N = f.shape[0]
    fe = np.concatenate([f,-f[::-1,:]],axis=0)
    df = np.real(ifft(ik*fft(fe,axis=0),axis=0))
    df = df[:N,:]
    
    return df

def ddx_none(f,ik):

    df = 0.
    
    return df

def SPECTRAL_x(sim):       # Set the differentiation operators
    
    if sim.Nx == 1:
        sim.ik = 0.
        
        sim.ddx_u = ddx_none
        sim.ddx_v = ddx_none
        sim.ddx_h = ddx_none
    elif sim.geomx == 'periodic':
        kx = 2*np.pi/sim.Lx*np.hstack([range(0,int(sim.Nx/2)), range(-int(sim.Nx/2),0)])
        sim.kx = kx.copy()
        sim.ik = 1j*np.tile(kx.reshape((sim.Nkx,1)),(1,sim.Ny))
        
        sim.ddx_u = ddx_period
        sim.ddx_v = ddx_period
        sim.ddx_h = ddx_period
    elif sim.geomx == 'walls':
        kx = np.pi/sim.Lx*np.hstack([range(0,int(sim.Nx)), range(-int(sim.Nx),0)])
        sim.kx = kx.copy()
        sim.ik = 1j*np.tile(kx.reshape((sim.Nkx,1)),(1,sim.Ny))
        
        sim.ddx_u = ddx_odd
        sim.ddx_v = ddx_even
        sim.ddx_h = ddx_even
    else:
        print "y boundary conditions must be from the list: periodic, walls"
        sys.exit()

        
