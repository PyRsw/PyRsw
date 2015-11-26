from __future__ import division
import numpy as np
import scipy.linalg as spalg
import matplotlib.pyplot as plt
import sys


## CHEB computes the Chebyshev differentiation matrix
## ------
#    matrix on N+1 points (i.e. N intervals)
#    D = differentiation matrix
#    x = Chebyshev grid

def cheb(N):
    if N == 0:
        D = 0
        x = 1
    else:
        x = np.cos(np.pi*np.array(range(0,N+1))/N).reshape([N+1,1])
        c = np.ravel(np.vstack([2, np.ones([N-1,1]), 2])) \
            *(-1)**np.ravel(np.array(range(0,N+1)))
        c = c.reshape(c.shape[0],1)
        X = np.tile(x,(1,N+1))
        dX = X-(X.conj().transpose())
        D  = (c*(1/c).conj().transpose())/(dX+(np.eye(N+1)))   # off-diagonal entries
        D  = D - np.diag(np.sum(D,1))   # diagonal entries
    return D,x
## ------



def stability_sw(sim, Dy, y, kk, ne):

    # Define matrices to store solution
    Nk = len(kk)
    growsw = np.zeros((ne,Nk))
    freqsw = np.zeros((ne,Nk))
    
    # Build matrices
    Z1 = np.zeros((sim.Ny+1,sim.Ny+1))
    I1 = np.eye(sim.Ny+1)

    # Set counter and start loop
    cnt = 0

    for k in kk:

        k2 = k**2

        # FJP: build part and then add k part for efficiency?
        # Set up Generalized Eigenvalue Problem (GEP)
        A = np.vstack((np.hstack((       np.diag(sim.UB,0),    (np.diag(sim.dUB - sim.f0))[:,1:-1],  sim.g*I1)),
                       np.hstack(( -sim.f0/k**2*I1[1:-1,:],           np.diag(sim.UB,0)[1:-1,1:-1], -sim.g/k**2*Dy[1:-1,:])),
                       np.hstack((         np.diag(sim.HB), (np.dot(Dy,np.diag(sim.HB,0)))[:,1:-1],  np.diag(sim.UB,0)))))

        # Solve Eigenvalue Problem Directly
        eigVals,eigVecs = spalg.eig(A)

        # Sort eigenvalues and eigenvectors
        ind = (-np.imag(eigVals)).argsort()
        eigVecs = eigVecs[ind,:]
        eigVals = k*eigVals[ind]

        growsw[:,cnt] = eigVals[0:ne].imag
        freqsw[:,cnt] = eigVals[0:ne].real

        cnt += 1

    return growsw, freqsw

