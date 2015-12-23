import numpy as np
import sys
 
def ddx_none(f,dx):
    return 0.

def ddx(f,dx):

    df = (f[1:,:] - f[0:-1,:])/dx
            
    return df

def ddx_periodic(f,dx):

    fs = np.concatenate([f[-1:,:],f,f[0:1,:]],axis=0)
    df = (fs[1:,:] - fs[0:-1,:])/dx
            
    return df

def ddx_walls(f,dx):

    fs = np.concatenate([-f[0:1,:],f,-f[-1:,:]],axis=0)
    df = (fs[1:,:] - fs[0:-1,:])/dx
            
    return df

def avx_none(f):
    return f


def avx(f):

    af = 0.5*(f[1:,:] + f[0:-1,:])
            
    return af

def avx_periodic(f):

    fs = np.concatenate([f[-1:,:],f,f[0:1,:]],axis=0)
    af = 0.5*(fs[1:,:] + fs[0:-1,:])
            
    return af

def avx_walls(f):

    fs = np.concatenate([-f[0:1,:],f,-f[-1:,:]],axis=0)
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

        try:
            import SADOURNY_x_fortran as fort

            sim.ddx_v = fort.ddx
            sim.ddx_h = fort.ddx
            sim.avx_v = fort.avx
            sim.avx_h = fort.avx

            if sim.geomx == 'periodic':
            
                sim.ddx_u = fort.ddx_periodic
                sim.avx_u = fort.avx_periodic

            elif sim.geomx == 'walls':

                sim.ddx_u = fort.ddx_walls
                sim.avx_u = fort.avx_walls

            else:
                print "x boundary conditions must be from the list: periodic, walls"
                sys.exit()

            print('Using Fortran in x')
        except:
            print('Unable to load Fortran for x functions')

            sim.ddx_v = ddx
            sim.ddx_h = ddx
            sim.avx_v = avx
            sim.avx_h = avx

            if sim.geomx == 'periodic':
            
                sim.ddx_u = ddx_periodic
                sim.avx_u = avx_periodic

            elif sim.geomx == 'walls':

                sim.ddx_u = ddx_walls
                sim.avx_u = avx_walls

            else:
                print "x boundary conditions must be from the list: periodic, walls"
                sys.exit()
