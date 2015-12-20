import numpy as np
import sys
 
def ddy_none(f,dy):

    df = 0.
    
    return df

def ddy(f,dy):

    df = (f[:,1:] - f[:,0:-1])/dy
            
    return df

def avy_none(f):

    af = f
            
    return af

def avy(f):

    af = 0.5*(f[:,1:] + f[:,0:-1])
            
    return af

def SADOURNY_y(sim):       # Set the differentiation operators

    if sim.Ny == 1:
        
        sim.ddy = ddy_none
        sim.avy = avy_none

    else:

        sim.ddy = ddy
        sim.avy = avy

