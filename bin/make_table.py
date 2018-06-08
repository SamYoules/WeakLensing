## SY 7/6/18
## Create a text file with a table of all unlensed and lensed RA & DEC for delta files

import fitsio
import numpy as N
import sys
import glob

# input directory name (lensed deltas)
dir = sys.argv[1]

# output file name
fout = open(sys.argv[2], 'w')

print('RA, DEC, Z, PLATE, MJD, FIBERID, THINGID, RA_LENS, DEC_LENS', file = fout)

alldeltas =  glob.glob(dir+'/*.fits.gz')
ndel = len(alldeltas)
i=0
for filename in alldeltas:
    hdus = fitsio.FITS(filename)
    print(i, ndel)
    i+=1

    for hdu in hdus[1:]:
        header = hdu.read_header()
        ra = header['RA']
        dec = header['DEC']
        ra_lens = header['RA_LENS']
        dec_lens = header['DEC_LENS']
        z = header['Z']
        plate = header['PLATE']
        mjd = header['MJD']
        fiberid = header['FIBERID']
        thingid = header['THING_ID']
        print (ra, dec, ra_lens, dec_lens, z, plate, mjd, fiberid, thingid, file = fout)

fout.close()
