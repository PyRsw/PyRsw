import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft, fftfreq, fftn, ifftn, fftshift

# Physical Parameters
day   = 86400.
lat   = 64.*(2.*np.pi/360)
Omega = 2*np.pi/day
f0    = 2.*Omega*np.sin(lat)
g0    = 0.5e-3
H0    = 4000.
c     = np.sqrt(g0*H0)
Rd    = c/(2.*Omega)

# Scale Parameters
alpha = 1.0
g = g0 / alpha
H = H0 * np.sqrt(alpha)

# Geometry Parameters
Nx, Ny = 64, 64
Lx, Ly = 400e3, 400e3
Lx, Ly = 2.*np.pi, 2.*np.pi
dx, dy = Lx/Nx, Ly/Ny

# Spatial Grid
x, y = np.arange(dx/2,Lx,dx), np.arange(dy/2,Ly,dy)
X, Y = np.meshgrid(x,y,indexing='ij')

# Spectral Grid
kx = 2*np.pi/Lx*np.hstack([range(0,int(Nx/2)), range(-int(Nx/2),0)])
ky = 2*np.pi/Ly*np.hstack([range(0,int(Ny/2)), range(-int(Ny/2),0)])
KX, KY = np.meshgrid(kx, ky, indexing='ij')
K = np.sqrt(KX**2 + KY**2)

print kx
print ky
print KX

K0 = kx[Nx/4]
W = 2*(kx[1] - kx[0])

theta = 2.*np.pi*np.random.randn(Nx,Ny)

etah = 1/K*np.exp(-((K-K0)/W)**2)*np.exp(1j*theta)
etah[0,0] = 0.0

eta = ifftn(etah).real

print K0
print W

plt.clf()
#plt.contourf(fftshift(KX),fftshift(KY), fftshift(np.abs(etah)))
plt.pcolormesh(X,Y,eta)
plt.colorbar()
plt.show()
