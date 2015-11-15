import numpy as np   
try:
    from pyfftw.interfaces.scipy_fftpack import fft, ifft, fftfreq
except:
    from scipy.fftpack import fft, ifft, fftfreq

def ddx_period(f,ik):

    df = np.real(ifft(ik*fft(f,axis=0),axis=0))

    return df

def ddx_even(f,ik):

    N = f.shape[0]
    fe = np.concatenate([f,f],axis=0)
    df = np.real(ifft(ik*fft(fe,axis=0),axis=0))
    df = df[:N,:]
    
    return df

def ddx_odd(f,ik):

    N = f.shape[0]
    fe = np.concatenate([f,-f],axis=0)
    df = np.real(ifft(ik*fft(fe,axis=0),axis=0))
    df = df[:N,:]
    
    return df


def SPECTRAL_x(sim):

    # Set the wavenumber vectors  
    if sim.geomx =='periodic':
        kx = 2*np.pi*fftfreq(sim.Nx,d=sim.Lx/sim.Nx)
        sim.kx = kx.copy()
        sim.ik = 1j*np.tile(kx.reshape((sim.Nx,1)),(1,sim.Ny))
    elif sim.geomx == 'walls':
        kx = np.pi*fftfreq(2*sim.Nx,d=sim.Lx/sim.Nx)
        sim.kx = kx.copy()
        sim.ik = 1j*np.tile(kx.reshape((2*sim.Nx,1)),(1,sim.Ny))
    else:
        print "x boundary conditions must be from the list: periodic, walls"
        sys.exit()

    # Set the differentiation operators
    if sim.geomx == 'periodic':
        sim.ddx_u = ddx_period
        sim.ddx_v = ddx_period
        sim.ddx_h = ddx_period
    elif sim.geomx == 'wall':
        sim.ddx_u = ddx_odd
        sim.ddx_v = ddx_even
        sim.ddx_h = ddx_even
