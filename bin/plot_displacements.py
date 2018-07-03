## SY 6/5/18
## Plots displacements of Mock QSOs when lensed by a giant blob placed in centre of NGC quasars

import numpy as N
import pylab as P
from astropy.table import Table
import sys
from collections import OrderedDict ## For removing duplicate labels in plot

# input directory name (e.g. outputs/table-blob.txt)
fin1 = sys.argv[1]

#P.rcParams.update({'font.size':22})

## Create a table of unlensed and lensed mocks qsos
t = Table(N.loadtxt(fin1, skiprows=1), names=('ra', 'dec', 'ra_lens', 'dec_lens', 'z', 'plate', 'mjd', 'fiberid', 'thing_id'))

## Pick random numbers to use as indices for plotting a thinned out sample of the data
rand_ind = N.random.choice(N.arange(len(t)), replace=False, size=1000)

P.figure(figsize=(8, 8))
for i in rand_ind:
    x = t['ra'][i]
    y = t['dec'][i]
    dx = t['ra_lens'][i] - t['ra'][i]
    dy = t['dec_lens'][i] - t['dec'][i]
    P.arrow(x, y, dx, dy, color='k', width=0.0003)
P.xlim((2.6,3.6))
P.ylim((0.0,1.0))
P.title('Displacements in Mock QSOs')
P.savefig('plots/sample_qsos.png')

