# -*- coding: utf-8 -*-
"""
Super basic script to solve heat diffusion equation using Euler's method
Use Dirichlet boundary condition with T = 0C at edge of grid

Created on Tue Feb 22 19:58:53 2022

@author: Gijs
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

print('starting...')

nx = 40        # number of cells
ny = 40        # number of cells
k = 384.1      # thermal conductivity of Cu in W·m−1·K−1
               # https://en.wikipedia.org/wiki/Thermal_conductivity
c = 0.385      # Cu specific heat capacity, J·g-1·K-1
rho = 8.96e6   # Cu density, g·m-3
dt = 0.2        # s
lx = 1         # m
ly = 1         # m
n_iter = 1000

all_x = np.linspace(0, lx, nx)
all_y = np.linspace(0, ly, ny)

cellsize_x = lx/nx
cellsize_y = ly/ny

xv, yv = np.meshgrid(all_x, all_y)

T0 = 200  # C
initial_shape = 'banana'
if initial_shape == 'boring_rectangle':
    # Initial temperature profile, just some hot square:
    T_init = T0 * ( ( xv >= 0.2 ) & ( xv <= 0.8 ) & ( yv >= 0.4 ) & ( yv <= 0.6 ) )
elif initial_shape == 'banana':
    # Alternative: a banana because why not
    T_init = T0 * ( ( ( (xv-0.5)**2 - yv ) < -0.2 ) & ( ( (xv-0.5)**2 - yv ) > -0.5 ) & ( yv >= 0.2 ) & ( yv <= 0.8 ) )

plt.figure()
plt.imshow(T_init, extent=[0,lx,0,ly])
plt.title('initial temperature profile')

# Run a single iteration
T = T_init.copy()
T_ = np.zeros(np.shape(xv))
for x in range(nx):
    for y in range(ny):
        if x == 0 or x == nx-1 or y == 0 or y == ny-1:
            T_[y, x] = 0
        else:
            T_[y, x] = T[y, x] - k*dt/(cellsize_y**2*c*rho) * (2 * T[y, x] - T[y-1, x] - T[y+1, x]) \
                               - k*dt/(cellsize_x**2*c*rho) * (2 * T[y, x] - T[y, x-1] - T[y, x+1] )
T = T_
plt.figure()
plt.imshow(T_, extent=[0,lx,0,ly])
plt.colorbar()
plt.title('after 1 iteration')

plt.figure()
plt.imshow(T_init - T_, extent=[0,lx,0,ly])
plt.colorbar()
plt.title('delta')


# Now run all iterations:
T = [ T_init ]
for i in range(1, n_iter+1):
    if np.mod(i, np.floor(n_iter/10)) == 0:
        print('iteration {} / {}...'.format(i, n_iter))
    T.append(np.zeros(np.shape(xv)))
    for x in range(nx):
        for y in range(ny):
            if x == 0 or x == nx-1 or y == 0 or y == ny-1:
                T[i][y, x] = 0
            else:
                T[i][y, x] = T[i-1][y, x] - k*dt/(cellsize_y**2*c*rho) * (2 * T[i-1][y, x] - T[i-1][y-1, x] - T[i-1][y+1, x]) \
                                          - k*dt/(cellsize_x**2*c*rho) * (2 * T[i-1][y, x] - T[i-1][y, x-1] - T[i-1][y, x+1] )

if False:
    for i in (0, 1, 2, 5, 10, 50, 100, 200, 500, 1000, 2000, 5000):
        plt.figure()
        plt.imshow(T[i], vmin=0, vmax=T0, extent=[0,lx,0,ly])
        plt.colorbar()
        plt.title('Iteration {}'.format(i))




# Generate animation
fig = plt.figure()
#creating a subplot 
ax1 = fig.add_subplot(1,1,1)

steps_per_frame = 20

def animate(i):
    ax1.clear()
    ax1.imshow(T[i*steps_per_frame], vmin=0, vmax=T0)
    plt.title('t = {:.1f} s'.format(i*dt*steps_per_frame))
    #plt.colorbar()
    
anim = animation.FuncAnimation(fig, animate, frames=int(n_iter/steps_per_frame), interval=20) 
anim.save('heat.gif', writer='imagemagick')
