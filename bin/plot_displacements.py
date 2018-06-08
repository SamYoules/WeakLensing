## SY 6/5/18
## Plots displacements of Mock QSOs when lensed by a giant blob placed in centre of NGC quasars

import numpy as N
import pylab as P
from astropy.table import Table

## Create a table of unlensed and lensed mocks qsos
t = Table(N.loadtxt('outputs/table-blob.txt', skiprows=1), names=('ra', 'dec', 'ra_lens', 'dec_lens', 'z', 'plate', 'mjd', 'fiberid', 'thing_id'))

P.figure()

P.plot(t['ra'], t['dec'], '.', color='k')
P.plot(t['ra_lens'], t['dec_lens'], 'r', color='r')  ## 
P.savefig('plots/all_qsos.png')

## Pick some random numbers to use as indices for plotting a thinned out sample of the data
rand_ind = N.random.choice(N.arange(len(t)), replace=False, size=1000)

P.figure(figsize=(8, 8))
for i in rand_ind:
    P.plot([t['ra'][i], t['ra_lens'][i]], [t['dec'][i], t['dec_lens'][i]], color='k')
    P.plot(t['ra_lens'][i], t['dec_lens'][i], '.', color='r')
P.xlim((2.6,3.6))
P.ylim((0.0,1.0))
#P.xlim((2,4))
#P.ylim((0.2,0.8))
P.title('Displacements in Mock QSOs')
P.savefig('plots/sample_qsos.png')
print(sum(t['ra_lens']-t['ra']!=0))
print(sum(t['ra_lens']-t['ra']==0))


