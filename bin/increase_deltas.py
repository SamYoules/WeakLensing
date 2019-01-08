# Author: Sam Youles
# Double the number of forests in each delta file, with RA & dec in same
# Healpixel, and rewrite to another directory

import numpy as N
import pylab as P
import scipy as sp
import fitsio
import healpy as hp
import glob
import sys
import os
from picca import constants
from picca.data import delta

indir = sys.argv[1]
outdir = sys.argv[2]

alldeltas = glob.glob(indir+"/*.fits.gz")
ndel = len(alldeltas)
i=1
for filename in alldeltas:
    hdus = fitsio.FITS(filename)
    print(i, ndel)
    i+=1

    fin = filename[-11:-8]
    out = fitsio.FITS(outdir + "/delta-" + fin + "-1.fits.gz",'rw',clobber=True)

    for d in hdus[1:]:
        header = d.read_header()
        ra = header['RA']
        dec = header['DEC']

        #-- Re-create columns 
        ll = d['LOGLAM'][:]
        de = d['DELTA'][:]
        we = d['WEIGHT'][:]
        co = d['CONT'][:] 
        cols=[ll, de, we, co]
        names=['LOGLAM','DELTA','WEIGHT','CONT']

        # make new random ra & dec within limits of pixel
        s = N.random.uniform(-1,1,2)*0.013

        #-- Rewrite new delta file with new values ### need to add forests OR create additional new delta
        header['RA'] = ra + s[0]
        header['DEC'] = dec + s[1]
        out.write(cols, names=names, header=header, \
              extname=str(header['THING_ID']))

    out.close()

