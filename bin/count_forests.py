import glob
import numpy as np
import healpy as hp
import fitsio

indir='deltas'
alldeltas = glob.glob(indir+'/*.fits.gz')
ndel = len(alldeltas)
i=0
j=0
for filename in alldeltas:
    hdus = fitsio.FITS(filename)
    i += len(hdus[1:])
    j += 1
    print(j, ndel)
print('number of forests = ', i)


