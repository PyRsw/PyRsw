##
## This file contains the core PyRsw user class.
##

import numpy as np
import matplotlib.pyplot as plt
import Plot_tools
import Diagnose
from scipy.fftpack import fftn, ifftn, fftfreq
import os, sys, shutil

def null_topo(a_sim):
    return

class Solution():
    def __init__(self,Nx,Ny,Nz):
        self.u = np.zeros((Nx,Ny,Nz+1))
        self.v = np.zeros((Nx,Ny,Nz+1))
        self.h = np.zeros((Nx,Ny,Nz+1))

class Flux():
    def __init__(self):
        self.u = []
        self.v = []
        self.h = []

class Simulation:

    # First-level initialization, default values
    def __init__(self):

        self.Nx = 1                 # Default grid size
        self.Ny = 1
        self.Nz = 1

        self.Lx = 1                 # Default domain sizes
        self.Ly = 1

        self.nfluxes = 0            # Length of flux history
        self.fluxes = Flux()
        
        self.method = 'Spectral'    # Spectral or Finite Volume (FV)?

        self.dynamics = 'Nonlinear' # Nonlinear or Linear

        self.num_threads = 1        # Number of threads for FFTW

        self.plot_vars = ['u','v','h']         # Which variables to plot
        self.ylims = []
        self.clims = []
        
        self.g    = 9.81            # gravity
        self.f0   = 1e-4            # Coriolis
        self.beta = 0.
        self.cfl  = 0.5             # default CFL
        self.time = 0               # initial time
        self.min_dt = 1e-3          # minimum timestep

        self.geomx = 'periodic'     # x boundary condition
        self.geomy = 'periodic'     # y boundary condition

        self.cmap = 'seismic'       # Default colour map
        
        self.run_name = 'test'      # Name of variable

        self.fcut = 0.6             # Filter cutoff
        self.ford = 2.0             # Filter order
        self.fstr = 20.0            # Filter strength

        self.vanishing = False
        self.fps = 15
        self.dpi = 150
        self.frame_count = 0
        self.out_counter = 0

        self.plott = np.inf
        self.diagt = np.inf

        self.num_steps = 0
        self.mean_dt = 0.

        self.topo_func = null_topo  # Default to no topograpy
        
    # Full initialization for once the user has specified parameters
    def initialize(self):

        print('Parameters:')
        print('-----------')
        if self.Nx > 1:
            print('geomx    = {0:s}'.format(self.geomx))
        if self.Ny > 1:
            print('geomy    = {0:s}'.format(self.geomy))
        print('stepper  = {0:s}'.format(self.stepper.__name__))
        print('method   = {0:s}'.format(self.method))
        print('dynamics = {0:s}'.format(self.dynamics))
        print('Nx       = {0:d}'.format(self.Nx))
        print('Ny       = {0:d}'.format(self.Ny))
        print('Nz       = {0:d}'.format(self.Nz))
        if self.f0 != 0:
            if self.beta != 0:
                print('Coriolis = beta-plane')
            else:
                print('Coriolis = f-plane')
        print ' '
        
        # Initialize grids and cell centres
        dxs = [1,1]
        
        if self.Nx > 1:
            dx = self.Lx/self.Nx
            self.x = np.arange(dx/2,self.Lx,dx) - self.Lx/2.
            dxs[0] = dx
        else:
            self.x = np.array([0.])
            
        if self.Ny > 1:
            dy = self.Ly/self.Ny
            self.y = np.arange(dy/2,self.Ly,dy) - self.Ly/2.
            dxs[1] = dy
        else:
            self.y = np.array([0.])
            
        self.dx = dxs

        if self.beta != 0.:
            if self.geomy != 'walls':
                print('beta-plane requires "walls" geometry in y.')
                sys.exit()

        # Define a 2D grid
        self.X, self.Y = np.meshgrid(self.x,self.y,indexing='ij')
        self.F = self.f0 + self.beta*self.Y

        # Initialize differentiation and averaging operators
        self.flux_method(self)
        self.x_derivs(self)
        self.y_derivs(self)

        # Compute reduced gravities
        if self.Nz > 1:
            tmp = np.ravel(self.rho.copy())
            tmp[1:] = tmp[1:] - tmp[0:-1]
            tmp = np.tile(tmp,(1,1)).reshape((1,len(self.rho)))
            self.gs = self.g*tmp/tmp[0,0]
        else:
            self.gs = np.array([[self.g]])

        # Initial conditions and topography
        self.soln = Solution(self.Nx,self.Ny,self.Nz)
        self.curr_flux = Solution(self.Nx,self.Ny,self.Nz)
        self.topo_func(self)

        # Construct the spectral filter
        fcut = self.fcut
        ford = self.ford
        fstr = self.fstr
        if self.Nx>1:
            k = self.kx/max(self.kx.ravel())
            filtx = np.exp(-fstr*((np.abs(k)-fcut)/(1-fcut))**ford)*(np.abs(k)>fcut) + (np.abs(k)<fcut)
            filtx = filtx.reshape((self.Nkx,1))
        else:
            filtx = np.array([1.0])
                
        if self.Ny>1:
            k = self.ky/max(self.ky.ravel())
            filty = np.exp(-fstr*((np.abs(k)-fcut)/(1-fcut))**ford)*(np.abs(k)>fcut) + (np.abs(k)<fcut)
            filty = filty.reshape((1,self.Nky))
        else:
            filty = np.array([1.0])
                
        self.sfilt = np.tile(filtx,(1,self.Nky))*np.tile(filty,(self.Nkx,1))
            
        
    def prepare_for_run(self):

        # If we're going to be diagnosing, initialize those
        Diagnose.initialize_diagnostics(self)
    
        # If we're saving, initialize those too
        if self.output:
            self.next_save_time = self.savet
    
        # If we're saving, initialize the directory
        if self.output or (self.animate == 'Save') or self.diagnose:
            self.initialize_saving()

        # If we're going to be plotting, then initialize the plots
        if self.animate != 'None':
            if (self.Nx > 1) and (self.Ny > 1):
                self.clims += [[]]*(len(self.plot_vars) - len(self.clims))
                self.initialize_plots = Plot_tools.initialize_plots_animsave_2D
            else:
                self.ylims += [[]]*(len(self.plot_vars) - len(self.ylims))
                self.initialize_plots = Plot_tools.initialize_plots_animsave_1D
        
            num_plot = self.end_time/self.plott+1
            if (self.Nx > 1) and (self.Ny == 1):
                self.hov_h = np.zeros((self.Nx,self.Nz,num_plot))
                self.hov_h[:,:,0] = self.soln.h[:,0,:-1]
            elif (self.Nx == 1) and (self.Ny > 1):
                self.hov_h = np.zeros((self.Ny,self.Nz,num_plot))
                self.hov_h[:,:,0] = self.soln.h[0,:,:-1]
            self.hov_count = 1

            self.initialize_plots(self)
            self.update_plots(self)
            self.next_plot_time = self.plott

    # Compute the current flux
    def flux(self):
        return self.flux_function(self)

    # Adjust time-step when necessary
    def adjust_dt(self):
        
        t = self.time + self.dt

        nt = self.end_time # next time for doing stuff
        if self.animate != 'None':
            nt = min([self.next_plot_time, nt])
        if self.diagnose:
            nt = min([self.next_diag_time, nt])
        if self.output:
            nt = min([self.next_save_time, nt])

        if nt < t:
            self.dt = nt - self.time

        t = self.time + self.dt

        do_plot = False
        if self.animate != 'None':
            if t + self.min_dt >= self.next_plot_time:
                do_plot = True

        do_diag = False
        if self.diagnose:
            if t + self.min_dt >= self.next_diag_time:
                do_diag = True

        do_save = False
        if self.output:
            if t + self.min_dt >= self.next_save_time:
                do_save = True

        return do_plot, do_diag, do_save

    # Advance the simulation one time-step.
    def step(self):

        self.compute_dt() 

        # Check if we need to adjust the time-step
        # to match an output time
        do_plot, do_diag, do_save = self.adjust_dt()

        self.stepper(self)

        # Filter if necessary
        if self.method == 'Spectral':
            self.apply_filter(self)

        self.time += self.dt
       
        if do_plot:
            self.update_plots(self)
            self.next_plot_time += self.plott
            if (self.Nx == 1) or (self.Ny == 1):
                Plot_tools.update_hov(self)

        if do_diag:
            Diagnose.update(self)

        if do_save:
            self.save_state()
            self.next_save_time += self.savet

        # Update the records
        self.mean_dt = (self.mean_dt*self.num_steps + self.dt)/(self.num_steps+1)
        self.num_steps += 1

        if do_plot or self.time == self.dt:

            h0 = self.soln.h[:,:,0] - self.soln.h[:,:,1]
            minh = np.min(np.ravel(h0))
            maxh = np.max(np.ravel(h0))

            maxu = np.max(np.ravel(self.soln.u[:,:,0]))
            minu = np.min(np.ravel(self.soln.u[:,:,0]))
            maxv = np.max(np.ravel(self.soln.v[:,:,0]))
            minv = np.min(np.ravel(self.soln.v[:,:,0]))
            
            mass = Diagnose.compute_mass(self)
            enrg = Diagnose.compute_PE(self) + Diagnose.compute_KE(self)

            if self.Nz == 1:
                pstr  = '{0:s}'.format(Plot_tools.smart_time(self.time))
                pstr += ' '*(13-len(pstr))
                try:
                    pstr += ',  dt = {0:0<7.1e}'.format(mean(self.dts))
                except:
                    pstr += ',  dt = {0:0<7.1e}'.format(self.dt)
                L = len(pstr) - 13
                pstr += ', min(u,v,h) = ({0: < 8.4e},{1: < 8.4e},{2: < 8.4e})'.format(minu,minv,minh)
                pstr += ', del_mass = {0: .2g}'.format(mass/self.Ms[0]-1)
                pstr += '\n'
                tmp = '  = {0:.3%}'.format(self.time/self.end_time)
                pstr += tmp
                pstr += ' '*(L - len(tmp))            
                pstr += 'avg = {0:0<7.1e}'.format(self.mean_dt)
                pstr += ', max(u,v,h) = ({0: < 8.4e},{1: < 8.4e},{2: < 8.4e})'.format(maxu,maxv,maxh)
                pstr += ', del_enrg = {0: .2g}'.format(enrg/(self.KEs[0]+self.PEs[0])-1)
            else:
                pstr  = 't = {0:.4g}s'.format(self.time)
                try:
                    pstr += ', dt = {0:0<7.1e}'.format(np.mean(self.dts))
                except:
                    pstr += ', dt = {0:0<7.1e}'.format(self.dt)
                pstr += ', del_mass = {0:+.2g}'.format(mass/self.Ms[0]-1)
                pstr += ', del_enrg = {0:+.2g}'.format(enrg/(self.KEs[0]+self.PEs[0])-1)
                pstr += '\n'
                pstr += '  = {0:.3%}'.format(self.time/self.end_time)

            head_str = ('\n{0:s}'.format(self.run_name))
            if self.animate != 'None':
                head_str += ': frame {0:d}'.format(self.frame_count-1)
            if self.output:
                head_str += ': output {0:d}'.format(self.out_counter-1)
            head_str += ': step {0:d}'.format(self.num_steps)

            print(head_str)
            print(pstr)

    # Step until end-time achieved.
    def run(self):

        self.prepare_for_run()

        while self.time < self.end_time:
            self.step()

        if self.diagnose:
            Diagnose.plot(self)
        
        if (self.animate == 'Anim'):
            plt.ioff()
            plt.show()

        if self.diagnose:
            Diagnose.save(self)

    # Compute time-step using CFL condition.
    def compute_dt(self):
        c = np.sqrt(self.g*np.sum(self.Hs))
        eps = 1e-8

        max_u = np.max(np.abs(self.soln.u.ravel()))
        max_v = np.max(np.abs(self.soln.v.ravel()))

        dt_x = self.end_time - ((self.Nx-1)/(self.Nx-1+eps))*(self.end_time - self.dx[0]/(max_u+2*c))
        dt_y = self.end_time - ((self.Ny-1)/(self.Ny-1+eps))*(self.end_time - self.dx[1]/(max_v+2*c))

        self.dt = max([self.cfl*min([dt_x,dt_y]),self.min_dt])

    # Initialize the saving
    def initialize_saving(self):
        
        path = 'Outputs/{0:s}'.format(self.run_name)
        if not(os.path.isdir('Outputs')):
            os.mkdir('Outputs')

        # If directory already exists, delete it.
        if os.path.isdir(path):
           print('Output directory {0:s} already exists. '.format(path) + \
                 'Warning, deleting everything in the directory.') 
           shutil.rmtree(path)
        # Make directory.
        os.mkdir(path)

        if self.animate == 'Save':
            os.mkdir(path+'/Frames')

        # Initialize saving stuff
        self.save_info()
        if self.output:
            self.save_state()
            self.save_grid()

    # Save information
    def save_info(self):
        fname = 'Outputs/{0:s}/info'.format(self.run_name)

        fp = open(fname, 'w')

        fp.write('Hs         = {0:s}\n'.format(self.array2str(self.Hs)))
        fp.write('rho        = {0:s}\n'.format(self.array2str(self.rho)))
        fp.write('g0         = {0:g}\n'.format(self.g))
        fp.write('Nx         = {0:d}\n'.format(self.Nx))
        fp.write('Ny         = {0:d}\n'.format(self.Ny))
        fp.write('Nz         = {0:d}\n'.format(self.Nz))
        if self.f0 > 0:
            fp.write('f0         = {0:g}\n'.format(self.f0))
        if self.beta > 0:
            fp.write('beta       = {0:g}\n'.format(self.beta))

        fp.close()

    # Save grid 
    def save_grid(self):
        fname = 'Outputs/{0:s}/grid'.format(self.run_name) 
        if self.Nx > 1 and self.Ny > 1:
            np.savez_compressed(fname, x = self.x, y = self.y, X = self.X, Y = self.Y)
        elif self.Nx > 1:
            np.savez_compressed(fname, x = self.x, X = self.X)
        elif self.Ny > 1:
            np.savez_compressed(fname, y = self.y, Y = self.Y)

    # Save current state
    def save_state(self):
        np.savez_compressed('Outputs/{0:s}/{1:05d}'.format(self.run_name,self.out_counter),
                    u = self.soln.u, v = self.soln.v, h = self.soln.h, t = self.time)

        self.out_counter += 1

    # Helper to write arrays
    def array2str(fp,arr):
        s = '[' + reduce(lambda s1,s2: s1+', '+s2, map(str,arr)) + ']'
        return s
