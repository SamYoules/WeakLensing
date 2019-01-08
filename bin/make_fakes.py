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

alldeltas = glob.glob("test_deltas/*.fits.gz")
ndel = len(alldeltas)
i=1
for filename in alldeltas:
    hdus = fitsio.FITS(filename)
    print(i, ndel)
    i+=1

    fin = filename[:-8]
    for j in N.arange(2,6):
        out = fitsio.FITS(fin + "-" + str(j) + ".fits.gz",'rw',clobber=True)

        for d in hdus[1:]:
            header = d.read_header()
            xi, yi, zi = sample_spherical(1)
            vect = N.asarray([xi,yi,zi])
            theta,phi = hp.vec2ang(vect)

            #-- Rewrite new delta file with new values
            header['RA'] = 2*sp.pi - phi[0]
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


