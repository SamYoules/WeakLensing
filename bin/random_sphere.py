# Plots an even distribution of points on sphere

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d

def sample_spherical(npoints, ndim=3):
    '''Creates 3D unit vector for a cube, then normalises to radius of sphere'''
    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)
    return vec


#-- Plot a grid on the sphere
phi = np.linspace(0, np.pi, 20)
theta = np.linspace(0, 2 * np.pi, 40)
x = np.outer(np.sin(theta), np.cos(phi))
y = np.outer(np.sin(theta), np.sin(phi))
z = np.outer(np.cos(theta), np.ones_like(phi))

#-- Get arrays of x,y,z coordinates for random points
xi, yi, zi = sample_spherical(500)

fig, ax = plt.subplots(1, 1, subplot_kw={'projection':'3d', 'aspect':'equal'})
ax.plot_wireframe(x, y, z, color='k', rstride=1, cstride=1)
ax.scatter(xi, yi, zi, s=100, c='r', zorder=10)

