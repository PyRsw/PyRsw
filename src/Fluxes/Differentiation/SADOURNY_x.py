import numpy as np
import sys
 
def ddx_none(f,dx):

    df = 0.0
    
    return df

def ddx(f,dx):

    df = (f[1:,:] - f[0:-1,:])/dx
            
    return df

def ddx_bdry(f,dx):

    fs = np.concatenate([f[-1:,:],f,f[0:1,:]],axis=0)
    df = (fs[1:,:] - fs[0:-1,:])/dx
            
    return df

def avx_none(f):

    af = f
            
    return af


def avx(f):

    af = 0.5*(f[1:,:] + f[0:-1,:])
            
    return af

def avx_bdry(f):

    fs = np.concatenate([f[-1:,:],f,f[0:1,:]],axis=0)
    af = 0.5*(fs[1:,:] + fs[0:-1,:])
            
    return af

def SADOURNY_x(sim):       # Set the differentiation operators

    if sim.Nx == 1:
        
        sim.ddx_u = ddx_none
        sim.ddx_v = ddx_none
        sim.ddx_h = ddx_none
        sim.avx_u = avx_none
        sim.avx_v = avx_none
        sim.avx_h = avx_none

    else:

        sim.ddx_u = ddx_bdry
        sim.ddx_v = ddx
        sim.ddx_h = ddx
        sim.avx_u = avx_bdry
        sim.avx_v = avx
        sim.avx_h = avx


        
