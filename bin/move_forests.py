# move_forests.py
# Author: Sam Youles
# 19/4/18
# Uses amap.fits (created in lensing_maps.py) to add bend angles from a lens to
# the RA and DEC of forests in the delta-{}.fits.gz files (from do_deltas.py)
# Writes textfile of old and new RA and DEC
# user input: run move_forests.py delta-files-directory-name

#from astropy.io import fits
import numpy as N
import healpy as hp
import sys
import glob
import os
import fitsio

# input
indir = sys.argv[1]
amap = sys.argv[2]
outdir = sys.argv[3]

# Load bend angles from alpha map into an array
amap_plus = hp.fitsfunc.read_map(amap, field=0)
amap_minus = hp.fitsfunc.read_map(amap, field=1)
NSIDE = 256

# Amend DEC and RA in each of the delta files by the bend angle from alpha map
alldeltas = glob.glob(indir+'/*.fits.gz')
ndel = len(alldeltas)
i=0
for filename in alldeltas:
    #hdus = fits.open(filename)
    hdus = fitsio.FITS(filename)
    print(i, ndel)
    i+=1

    out = fitsio.FITS(outdir+"/"+os.path.basename(filename),'rw',clobber=True)

    for hdu in hdus[1:]:
        header = hdu.read_header()
        ra = header['RA']
        dec = header['DEC']

        # Convert RA and DEC to index to find corresponding pixel in amap
        # Note: theta is defined in the range [0,pi] and declination is defined in
        # the range [-pi/2,pi/2].
        index = hp.pixelfunc.ang2pix(NSIDE, -dec + N.pi/2., ra)

        # Add bend angles to ra and dec
        ra_lens = ra - amap_plus[index]
        dec_lens = dec - amap_minus[index]

        # Rewrite new delta file with new values
        header['RA_LENS'] = ra_lens
        header['DEC_LENS'] = dec_lens
      
        #-- Re-create columns (maybe there's a better way to do this?) 
        ll = hdu['LOGLAM'][:]
        de = hdu['DELTA'][:]
        we = hdu['WEIGHT'][:]
        co = hdu['CONT'][:] 
        cols=[ll, de, we, co]
        names=['LOGLAM','DELTA','WEIGHT','CONT']
        out.write(cols, names=names, header=header, \
                  extname=str(header['THING_ID']))

    out.close()

