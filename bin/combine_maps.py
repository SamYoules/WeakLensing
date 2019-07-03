#!/usr/bin/env python

##- SY 4/2/19
##- Make kappa map from weighted average of 2 maps
##- sys args are:
##-    auto map name
##-    cross map name
##-    output map name

import numpy as N
import pylab as P
import scipy as sp
import fitsio
from astropy.io import fits
import healpy
import sys

##- User input (auto and cross kappa maps)
map1   = sys.argv[1]
map2   = sys.argv[2]
mapout = sys.argv[3]

##- Read maps (kappa and weight)
khead = fits.open(map1)[1].header
k1 = fits.open(map1)[1].data.kappa
w1 = fits.open(map1)[1].data.wkappa
k2 = fits.open(map2)[1].data.kappa
w2 = fits.open(map2)[1].data.wkappa

##- Calculate weighted average
kap = N.zeros(w1.size)
wkap = N.zeros(w1.size)

for i in range(w1.size):
    if (w1[i]==0):
        kap[i] = 0.
        wkap[i] = 0.
    else:
        kap[i] = (k1[i]*w1[i] + k2[i]*w2[i])/(w1[i] + w2[i])
        wkap[i] = w1[i] + w2[i]

##- Write map
out = fitsio.FITS(mapout, 'rw', clobber=True)
head = {}
head['RPMIN']=khead['RPMIN']
head['RPMAX']=khead['RPMAX']
head['RTMIN']=khead['RTMIN']
head['RTMAX']=khead['RTMAX']
head['NT']=khead['NT']
head['NP']=khead['NP']
head['NSIDE']=khead['NSIDE']
out.write([kap, wkap], names=['kappa', 'wkappa'], header=head)
out.close()

