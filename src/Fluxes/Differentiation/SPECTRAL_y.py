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
    fe = np.concatenate([f,f[:,::-1]],axis=1)
    fe = np.concatenate([f,f],axis=1)
    df = np.real(ifft(il*fft(fe,axis=1),axis=1))
    df = df[:,:N]
    
    return df

def ddy_odd(f,il):

    N = f.shape[1]
    fe = np.concatenate([f,-f[:,::-1]],axis=1)
    df = np.real(ifft(il*fft(fe,axis=1),axis=1))
    df = df[:,:N]
    
    return df

def ddy_none(f,ik):

    df = 0.
    
    return df

def SPECTRAL_y(sim):  # Set the differentiation operators
    
    if sim.Ny == 1:
        sim.il = 0
        
        sim.ddy_u = ddy_none
        sim.ddy_v = ddy_none
        sim.ddy_h = ddy_none
    elif sim.geomy == 'periodic':
        ky = 2*np.pi/sim.Ly*np.hstack([range(0,int(sim.Ny/2)), range(-int(sim.Ny/2),0)])
        sim.ky = ky.copy()
        sim.il = 1j*np.tile(ky.reshape((1,sim.Nky)),(sim.Nkx,1))
        
        sim.ddy_u = ddy_period
        sim.ddy_v = ddy_period
        sim.ddy_h = ddy_period
    elif sim.geomy == 'walls':
        ky = np.pi/sim.Ly*np.hstack([range(0,int(sim.Ny)), range(-int(sim.Ny),0)])
        sim.ky = ky.copy()
        sim.il = 1j*np.tile(ky.reshape((1,sim.Nky)),(sim.Nkx,1))
        
        sim.ddy_u = ddy_even
        sim.ddy_v = ddy_odd
        sim.ddy_h = ddy_even
    else:
        print "y boundary conditions must be from the list: periodic, walls"
        sys.exit()
