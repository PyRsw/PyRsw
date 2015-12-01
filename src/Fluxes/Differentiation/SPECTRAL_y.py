import numpy as np   
try:
    import pyfftw
    from pyfftw.interfaces.scipy_fftpack import fft, ifft, fftfreq
except:
    from scipy.fftpack import fft, ifft, fftfreq


def ddy_none(f,ik):

    df = 0.
    
    return df


def SPECTRAL_y(sim):       # Set the differentiation operators

    if sim.Ny == 1:
        sim.il = 0.
        
        sim.ddy_u = ddy_none
        sim.ddy_v = ddy_none
        sim.ddy_h = ddy_none

    else:

        # If possible, use pyfftw to preserve wisdom.
        try:

            tmp_in_u  = pyfftw.n_byte_align_empty((sim.Nx,sim.Nky),16,dtype='complex128')
            tmp_out_u = pyfftw.n_byte_align_empty((sim.Nx,sim.Nky),16,dtype='complex128')
            sim.ffty_u  = pyfftw.FFTW(tmp_in_u, tmp_out_u, axes=[1], direction='FFTW_FORWARD', threads = sim.num_threads)
            sim.iffty_u = pyfftw.FFTW(tmp_out_u, tmp_in_u, axes=[1], direction='FFTW_BACKWARD', threads = sim.num_threads)

            tmp_in_v  = pyfftw.n_byte_align_empty((sim.Nx,sim.Nky),16,dtype='complex128')
            tmp_out_v = pyfftw.n_byte_align_empty((sim.Nx,sim.Nky),16,dtype='complex128')
            sim.ffty_v  = pyfftw.FFTW(tmp_in_v, tmp_out_v, axes=[1], direction='FFTW_FORWARD', threads  = sim.num_threads)
            sim.iffty_v = pyfftw.FFTW(tmp_out_v, tmp_in_v, axes=[1], direction='FFTW_BACKWARD', threads = sim.num_threads)

            tmp_in_h  = pyfftw.n_byte_align_empty((sim.Nx,sim.Nky),16,dtype='complex128')
            tmp_out_h = pyfftw.n_byte_align_empty((sim.Nx,sim.Nky),16,dtype='complex128')
            sim.ffty_h  = pyfftw.FFTW(tmp_in_h, tmp_out_h, axes=[1], direction='FFTW_FORWARD', threads  = sim.num_threads)
            sim.iffty_h = pyfftw.FFTW(tmp_out_h, tmp_in_h, axes=[1], direction='FFTW_BACKWARD', threads = sim.num_threads)

        except:
            def iffty(ar):
                return ifft(ar,axis=1)
            def ffty(ar):
                return fft(ar,axis=1)

            sim.ffty_u  = ffty
            sim.iffty_u = iffty
            sim.ffty_v  = ffty
            sim.iffty_v = iffty
            sim.ffty_h  = ffty
            sim.iffty_h = iffty
            print('Failed to use fftw for y transforms')

        if sim.geomy == 'periodic':

            def ddy_u(f,sim):
                df = np.real(sim.iffty_u(sim.il*sim.ffty_u(f)))
                return df
    
            def ddy_v(f,sim):
                df = np.real(sim.iffty_v(sim.il*sim.ffty_v(f)))
                return df
    
            def ddy_h(f,sim):
                df = np.real(sim.iffty_h(sim.il*sim.ffty_h(f)))
                return df

            ky = 2*np.pi/sim.Ly*np.hstack([range(0,int(sim.Ny/2)), range(-int(sim.Ny/2),0)])
            sim.ky = ky.copy()
            sim.il = 1j*np.tile(ky.reshape((1,sim.Nky)),(sim.Nx,1))
    
        elif sim.geomy == 'walls':

            def ddy_u(f,sim):
                N = f.shape[1]
                fe = np.concatenate([f,f[:,::-1]],axis=1)
                df = np.real(sim.iffty_u(sim.il*sim.ffty_u(fe)))[:,:N]
                return df
    
            def ddy_v(f,sim):
                N = f.shape[1]
                fe = np.concatenate([f,-f[:,::-1]],axis=1)
                df = np.real(sim.iffty_v(sim.il*sim.ffty_v(fe)))[:,:N]
                return df
    
            def ddy_h(f,sim):
                N = f.shape[1]
                fe = np.concatenate([f,f[:,::-1]],axis=1)
                df = np.real(sim.iffty_h(sim.il*sim.ffty_h(fe)))[:,:N]
                return df

            ky = np.pi/sim.Ly*np.hstack([range(0,int(sim.Ny)), range(-int(sim.Ny),0)])
            sim.ky = ky.copy()
            sim.il = 1j*np.tile(ky.reshape((1,sim.Nky)),(sim.Nx,1))
        
        else:
            print "y boundary conditions must be from the list: periodic, walls"
            sys.exit()

        sim.ddy_u = ddy_u
        sim.ddy_v = ddy_v
        sim.ddy_h = ddy_h


