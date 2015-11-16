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

def filter_periodic(sim):
    for ii in range(sim.Nz):

        sim.soln.u[:,:,ii] = ifftn(sim.sfilt*fftn(sim.soln.u[:,:,ii],axes=[0,1]),axes=[0,1]).real
        sim.soln.v[:,:,ii] = ifftn(sim.sfilt*fftn(sim.soln.v[:,:,ii],axes=[0,1]),axes=[0,1]).real
        sim.soln.h[:,:,ii] = ifftn(sim.sfilt*fftn(sim.soln.h[:,:,ii],axes=[0,1]),axes=[0,1]).real

        
def filter_walls(sim):
    for ii in range(sim.Nz):

        fe =  np.concatenate([sim.soln.u[:,:,ii],-sim.soln.u[:,:,ii]],axis=0)
        fe = ifftn(sim.sfilt*fftn(fe,axes=[0,1]),axes=[0,1]).real
        sim.soln.u[:,:,ii] = fe[:sim.Nx,:]

        fe =  np.concatenate([sim.soln.v[:,:,ii],sim.soln.v[:,:,ii]],axis=0)
        fe = ifftn(sim.sfilt*fftn(fe,axes=[0,1]),axes=[0,1]).real
        sim.soln.v[:,:,ii] = fe[:sim.Nx,:]
                
        fe =  np.concatenate([sim.soln.h[:,:,ii],sim.soln.h[:,:,ii]],axis=0)
        fe = ifftn(sim.sfilt*fftn(fe,axes=[0,1]),axes=[0,1]).real
        sim.soln.h[:,:,ii] = fe[:sim.Nx,:]

                
def SPECTRAL_x(sim):

    # Set the wavenumber vectors  
    if sim.geomx =='periodic':
        kx = 2*np.pi/sim.Lx*np.hstack([range(0,int(sim.Nx/2)), range(-int(sim.Nx/2),0)])
        sim.kx = kx.copy()
        sim.ik = 1j*np.tile(kx.reshape((sim.Nx,1)),(1,sim.Ny))
    elif sim.geomx == 'walls':
        kx = np.pi/sim.Lx*np.hstack([range(0,int(sim.Nx)), range(-int(sim.Nx),0)])
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
        sim.apply_filter = filter_periodic
    elif sim.geomx == 'walls':    #FJP: why does v satisfy Dirichlet BCs and not u?
        sim.ddx_u = ddx_odd
        sim.ddx_v = ddx_even
        sim.ddx_h = ddx_even
        sim.apply_filter = filter_walls
    else:
        print "x boundary conditions must be from the list: periodic, walls"
        sys.exit()
        
