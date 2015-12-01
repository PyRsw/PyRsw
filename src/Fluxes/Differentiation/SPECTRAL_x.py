import numpy as np
import sys
 
try:
    import pyfftw
    from pyfftw.interfaces.scipy_fftpack import fft, ifft, fftfreq, fftn, ifftn
except:
    from scipy.fftpack import fft, ifft, fftfreq, fftn, ifftn

def ddx_none(f,ik):

    df = 0.
    
    return df

def SPECTRAL_x(sim):       # Set the differentiation operators

    if sim.Nx == 1:
        sim.ik = 0.
        
        sim.ddx_u = ddx_none
        sim.ddx_v = ddx_none
        sim.ddx_h = ddx_none

    else:

        # If possible, use pyfftw to preserve wisdom.
        try:

            tmp_in_u  = pyfftw.n_byte_align_empty((sim.Nkx,sim.Ny),16,dtype='complex128')
            tmp_out_u = pyfftw.n_byte_align_empty((sim.Nkx,sim.Ny),16,dtype='complex128')
            sim.fftx_u  = pyfftw.FFTW(tmp_in_u, tmp_out_u, axes=[0], direction='FFTW_FORWARD', threads = sim.num_threads)
            sim.ifftx_u = pyfftw.FFTW(tmp_out_u, tmp_in_u, axes=[0], direction='FFTW_BACKWARD', threads = sim.num_threads)

            tmp_in_v  = pyfftw.n_byte_align_empty((sim.Nkx,sim.Ny),16,dtype='complex128')
            tmp_out_v = pyfftw.n_byte_align_empty((sim.Nkx,sim.Ny),16,dtype='complex128')
            sim.fftx_v  = pyfftw.FFTW(tmp_in_v, tmp_out_v, axes=[0], direction='FFTW_FORWARD', threads = sim.num_threads)
            sim.ifftx_v = pyfftw.FFTW(tmp_out_v, tmp_in_v, axes=[0], direction='FFTW_BACKWARD', threads = sim.num_threads)

            tmp_in_h  = pyfftw.n_byte_align_empty((sim.Nkx,sim.Ny),16,dtype='complex128')
            tmp_out_h = pyfftw.n_byte_align_empty((sim.Nkx,sim.Ny),16,dtype='complex128')
            sim.fftx_h  = pyfftw.FFTW(tmp_in_h, tmp_out_h, axes=[0], direction='FFTW_FORWARD', threads = sim.num_threads)
            sim.ifftx_h = pyfftw.FFTW(tmp_out_h, tmp_in_h, axes=[0], direction='FFTW_BACKWARD', threads = sim.num_threads)

        except:
            def ifftx(ar):
                return ifft(ar,axis=0)
            def fftx(ar):
                return fft(ar,axis=0)

            sim.fftx_u  = fftx
            sim.ifftx_u = ifftx
            sim.fftx_v  = fftx
            sim.ifftx_v = ifftx
            sim.fftx_h  = fftx
            sim.ifftx_h = ifftx
            print('Failed to use fftw for x transforms')

        if sim.geomx == 'periodic':

            def ddx_u(f,sim):
                df = np.real(sim.ifftx_u(sim.ik*sim.fftx_u(f)))
                return df
    
            def ddx_v(f,sim):
                df = np.real(sim.ifftx_v(sim.ik*sim.fftx_v(f)))
                return df
    
            def ddx_h(f,sim):
                df = np.real(sim.ifftx_h(sim.ik*sim.fftx_h(f)))
                return df
    
            kx = 2*np.pi/sim.Lx*np.hstack([range(0,int(sim.Nx/2)), range(-int(sim.Nx/2),0)])
            sim.kx = kx.copy()
            sim.ik = 1j*np.tile(kx.reshape((sim.Nkx,1)),(1,sim.Ny))
    
        elif sim.geomx == 'walls':

            def ddx_u(f,sim):
                N = f.shape[0]
                fe = np.concatenate([f,-f[::-1,:]],axis=0)
                df = np.real(sim.ifftx_u(sim.ik*sim.fftx_u(fe)))[:N,:]
                return df
    
            def ddx_v(f,sim):
                N = f.shape[0]
                fe = np.concatenate([f,f[::-1,:]],axis=0)
                df = np.real(sim.ifftx_v(sim.ik*sim.fftx_v(fe)))[:N,:]
                return df
    
            def ddx_h(f,sim):
                N = f.shape[0]
                fe = np.concatenate([f,f[::-1,:]],axis=0)
                df = np.real(sim.ifftx_h(sim.ik*sim.fftx_h(fe)))[:N,:]
                return df

            kx = np.pi/sim.Lx*np.hstack([range(0,int(sim.Nx)), range(-int(sim.Nx),0)])
            sim.kx = kx.copy()
            sim.ik = 1j*np.tile(kx.reshape((sim.Nkx,1)),(1,sim.Ny))
        
        else:
            print "x boundary conditions must be from the list: periodic, walls"
            sys.exit()

        sim.ddx_u = ddx_u
        sim.ddx_v = ddx_v
        sim.ddx_h = ddx_h


        
