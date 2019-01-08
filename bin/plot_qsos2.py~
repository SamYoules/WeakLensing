import numpy as N
from astropy.table import Table
import sys
from collections import OrderedDict ## For removing duplicate labels in plot
import healpy as hp
from matplotlib import pyplot as plt

# input directory name (e.g. wsky_table.txt)
fin1 = sys.argv[1]

## Create a table of unlensed and lensed mocks qsos
t = Table(N.loadtxt(fin1, skiprows=1), names=('ra', 'dec', 'ra_lens', 'dec_lens', 'z', 'plate', 'mjd', 'fiberid', 'thing_id'))

## Pick random numbers to use as indices for plotting a thinned out sample of the data
rand_ind = N.random.choice(N.arange(len(t)), replace=False, size=5000)

th=[]
ph=[]

for i in rand_ind:
    th.append(N.pi/2. - t['dec'][i])
    ph.append(t['ra'][i])

hp.mollview()
theta=N.asarray(th)
phi = N.asarray(ph)
hp.projscatter(theta,phi)
plt.show()


