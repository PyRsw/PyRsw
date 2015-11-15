import numpy as np
import WENOx_periodic
import matplotlib.pyplot as plt
Nx = 20
x = np.linspace(0,1,Nx+1)[:-1].reshape((Nx,1))
y = np.sin(2*np.pi*x)
dye = 2*np.pi*np.cos(2*np.pi*x)
dx = [x[1,0]-x[0,0],1]
dym, dyp = WENOx_periodic.ddx(y,dx)
xp = x + dx[0]/2.
xm = x - dx[0]/2.
plt.plot(xm,dym,'-r',x,dye,'-b',xp,dyp,'-g')
plt.show()
