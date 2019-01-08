# Author: Sam Youles
# Remove a randomised proportion of forests from each delta file and rewrite to
# another directory

import numpy as np
import healpy as hp
import sys
import glob
import os
import fitsio
from kappa_lya import *

# input directory name containing delta files
indir = sys.argv[1]
outdir = sys.argv[2]

seed=1
np.random.seed(seed)

# Remove a randomised proportion of forests from each delta file and rewrite to another directory
alldeltas = glob.glob(indir+'/*.fits.gz')
ndel = len(alldeltas)
i=0
for filename in alldeltas:
    hdus = fitsio.FITS(filename)
    print(i, ndel)
    i+=1

    out = fitsio.FITS(outdir+"/"+os.path.basename(filename),'rw',clobber=True)

    for hdu in hdus[1:]:
        header = hdu.read_header()
        ll = hdu['LOGLAM'][:]
        de = hdu['DELTA'][:]
        we = hdu['WEIGHT'][:]
        co = hdu['CONT'][:] 
        a = (len(ll))
        rand_numbers = np.random.rand(len(ll))
        w = (rand_numbers<0.27)
        ll = ll[w]
        de = de[w]
        we = we[w]
        co = co[w]
        print(a, len(ll))
        cols=[ll, de, we, co]
        names=['LOGLAM','DELTA','WEIGHT','CONT']
        out.write(cols, names=names, header=header, \
                  extname=str(header['THING_ID']))

    out.close()

