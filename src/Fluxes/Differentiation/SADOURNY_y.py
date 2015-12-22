import numpy as np
import sys
 
def ddy_none(f,dy):

    df = 0.0
    
    return df

def ddy(f,dy):

    df = (f[:,1:] - f[:,0:-1])/dy
            
    return df

def ddy_periodic(f,dy):

    fs = np.concatenate([f[:,-1:],f,f[:,0:1]],axis=1)
    df = (fs[:,1:] - fs[:,0:-1])/dy
            
    return df

def ddy_walls(f,dy):

    fs = np.concatenate([-f[:,0:1],f,-f[:,-1:]],axis=1)
    df = (fs[:,1:] - fs[:,0:-1])/dy
            
    return df

def avy_none(f):

    af = f
            
    return af

def avy(f):

    af = 0.5*(f[:,1:] + f[:,0:-1])
            
    return af

def avy_periodic(f):

    fs = np.concatenate([f[:,-1:],f,f[:,0:1]],axis=1)
    af = 0.5*(fs[:,1:] + fs[:,0:-1])
            
    return af

def avy_walls(f):

    fs = np.concatenate([-f[:,0:1],f,-f[:,-1:]],axis=1)
    af = 0.5*(fs[:,1:] + fs[:,0:-1])
    
    return af

def SADOURNY_y(sim):       # Set the differentiation operators

    if sim.Ny == 1:
        
        sim.ddy_u = ddy_none
        sim.ddy_v = ddy_none
        sim.ddy_h = ddy_none
        sim.avy_u = avy_none
        sim.avy_v = avy_none
        sim.avy_h = avy_none

    else:

        sim.ddy_u = ddy
        sim.ddy_h = ddy
        sim.avy_u = avy
        sim.avy_h = avy

        if sim.geomy == 'periodic':
            
            sim.ddy_v = ddy_periodic
            sim.avy_v = avy_periodic
            
        elif sim.geomy == 'walls':
            
            sim.ddy_v = ddy_walls
            sim.avy_v = avy_walls

        else:
            print "y boundary conditions must be from the list: periodic, walls"
            sys.exit()


            
