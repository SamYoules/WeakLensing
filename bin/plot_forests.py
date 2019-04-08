import glob
import fitsio
import pylab as P
import numpy as N

f2 = glob.glob("deltas_noiseless/*.fits.gz")
ndel = len(f2)
i=1
ra = []
dec = []

for filename in f2:
    print(i, ndel)
    i+=1
    hdus = fitsio.FITS(filename)
    for d in hdus[1:]:
        header = d.read_header()
        ra.append(header['RA']/N.pi * 180)
        dec.append(header['DEC']/N.pi * 180)

P.figure()
P.scatter(ra,dec)
