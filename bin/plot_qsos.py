import numpy as N
#import pylab as P
from astropy.table import Table
import sys
from collections import OrderedDict ## For removing duplicate labels in plot
import healpy as hp
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import healpy as hp

# input directory name (e.g. wsky_table.txt)
fin1 = sys.argv[1]

#P.rcParams.update({'font.size':20})

## Create a table of unlensed and lensed mocks qsos
t = Table(N.loadtxt(fin1, skiprows=1), names=('ra', 'dec', 'ra_lens', 'dec_lens', 'z', 'plate', 'mjd', 'fiberid', 'thing_id'))

## Pick random numbers to use as indices for plotting a thinned out sample of the data
rand_ind = N.random.choice(N.arange(len(t)), replace=False, size=1000)

#-- Plot a grid on the sphere
phi = N.linspace(0, N.pi, 20)
theta = N.linspace(0, 2 * N.pi, 40)
x = N.outer(N.sin(theta), N.cos(phi))
y = N.outer(N.sin(theta), N.sin(phi))
z = N.outer(N.cos(theta), N.ones_like(phi))

fig, ax = plt.subplots(1, 1, subplot_kw={'projection':'3d', 'aspect':'equal'})
ax.plot_wireframe(x, y, z, color='k', rstride=1, cstride=1)

#-- Plot qsos
for i in rand_ind:
    ph = t['ra'][i]
    th = N.pi/2. - t['dec'][i]
    x1,x2,x3=(hp.ang2vec(th,ph))
    ax.scatter(x1, x2, x3, s=50, c='r', zorder=10)
plt.show()

