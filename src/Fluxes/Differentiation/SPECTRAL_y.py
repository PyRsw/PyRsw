import numpy as np   
try:
    from pyfftw.interfaces.scipy_fftpack import fft, ifft, fftfreq
except:
    from scipy.fftpack import fft, ifft, fftfreq

def ddy_period(f,il):

    df = np.real(ifft(il*fft(f,axis=1),axis=1))

    return df

def ddy_even(f,il):

    N = f.shape[1]
    fe = np.concatenate([f,f],axis=1)
    df = np.real(ifft(il*fft(fe,axis=1),axis=1))
    df = df[:,:N]
    
    return df

def ddy_odd(f,il):

    N = f.shape[1]
    fe = np.concatenate([f,-f],axis=1)
    df = np.real(ifft(il*fft(fe,axis=1),axis=1))
    df = df[:,:N]
    
    return df


def SPECTRAL_y(sim):

    # Set the wavenumber vectors  
    if sim.geomy =='periodic':
        ky = 2*np.pi*fftfreq(sim.Ny,d=sim.Ly/sim.Ny)
        sim.ky = ky.copy()
        sim.il = 1j*np.tile(ky.reshape((1,sim.Ny)),(sim.Nx,1))
    elif sim.geomy == 'walls':
        ky = 1j*np.pi*fftfreq(2*sim.Ny,d=sim.Ly/sim.Ny)
        sim.ky = ky.copy()
        sim.il = 1j*np.tile(ky.reshape((1,2*sim.Ny)),(sim.Nx,1))
    else:
        print "x boundary conditions must be from the list: periodic, walls"
        sys.exit()

    # Set the differentiation operators
    if sim.geomy == 'periodic':
        sim.ddy_u = ddy_period
        sim.ddy_v = ddy_period
        sim.ddy_h = ddy_period
    elif sim.geomy == 'wall':
        sim.ddy_u = ddy_even
        sim.ddy_v = ddy_odd
        sim.ddy_h = ddy_even
