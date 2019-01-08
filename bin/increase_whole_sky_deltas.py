## SY 30/7/18
## Writes new delta files with randomised RA & dec but otherwise with same
## data as original (original folder contained deltas suffixed 1,2,3,4 or 5
## New files will have suffixes 6,7,8,9 or 0.

import numpy as N
import pylab as P
import scipy as sp
import fitsio
import healpy as hp
import glob
from picca import constants
from picca.data import delta

def sample_spherical(npoints, ndim=3):
    '''Get coordinates for a random point on a sphere'''
    vec = N.random.randn(ndim, npoints)
    vec /= N.linalg.norm(vec, axis=0)
    return vec

alldeltas = glob.glob("ws_deltas/*.fits.gz")
ndel = len(alldeltas)
i=1
for filename in alldeltas:
    hdus = fitsio.FITS(filename)
    print(i, ndel)
    i+=1

    fin = filename[:-9]
    k = filename[-9:-8]
    if k == '5':
        j = 0
    else:
        j = int(k) + 5
    out = fitsio.FITS(fin + str(j) + ".fits.gz",'rw',clobber=True)

    for d in hdus[1:]:
        header = d.read_header()
        xi, yi, zi = sample_spherical(1)
        vect = N.asarray([xi,yi,zi])
        theta,phi = hp.vec2ang(vect)

        #-- Rewrite new delta file with new values
        header['RA'] = phi[0]
        header['DEC'] = sp.pi/2. - theta[0]

        #-- Re-create columns 
        ll = d['LOGLAM'][:]
        de = d['DELTA'][:]
        we = d['WEIGHT'][:]
        co = d['CONT'][:] 
        cols=[ll, de, we, co]
        names=['LOGLAM','DELTA','WEIGHT','CONT']
        out.write(cols, names=names, header=header, \
              extname=str(header['THING_ID']))

    out.close()


