import numpy as np
import sys
 
def ddx_none(f,dx):

    df = 0.
    
    return df

def ddx(f,dx):

    df = (f[1:,:] - f[0:-1,:])/dx
            
    return df

def avx_none(f):

    af = f
            
    return af


def avx(f):

    af = 0.5*(f[1:,:] + f[0:-1,:])
            
    return af


def SADOURNY_x(sim):       # Set the differentiation operators

    if sim.Nx == 1:
        
        sim.ddx = ddx_none
        sim.avx = avx_none

    else:

        sim.ddx = ddx
        sim.avx = avx


        
